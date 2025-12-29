#!/usr/bin/env python3

from time import sleep
from subprocess import CalledProcessError, SubprocessError
from sys import exit
from pathlib import Path
import subprocess as subp
from random import randint
import multiprocessing as mp
from collections import deque
from shutil import which, move, rmtree
from typing import Union as OneOf, List
from argparse import ArgumentParser, Namespace

MAJOR_VERSION_NUMBER = 0
MINOR_VERSION_NUMBER = 1
PATCH_VERSION_NUMBER = 3

VERSION_STRING = f"v{MAJOR_VERSION_NUMBER}.{MINOR_VERSION_NUMBER}.{PATCH_VERSION_NUMBER}"

# script sets the 0x400 tag for optical and PCR duplicates
# (PCR duplicates should actually not occur due to UMI wetlab protocol
#  and in silico collapsing of this information)

rsync_path = which('rsync')
if not Path(rsync_path).is_file():
    raise FileNotFoundError(f"the 'rsync' executable path was not found or is not accessible")


def get_cmdline_args() -> Namespace:
    parser = ArgumentParser()
    # versioning
    parser.add_argument('--version', action='version',
                        version=f"{VERSION_STRING}")
    parser.add_argument('-tmp', '--temp-dir', metavar='File', dest="temp_dir", type=Path, required=True,
                        help='The temporary directory for the scaffold-wise correction. Scaffold BAM files '
                             'are concatenated afterwards to form the final BAM file - so twice the inputs '
                             'size is required for this operation. Temporary directory will be deleted after '
                             'rsyncing the result back to the target directory.')
    parser.add_argument('-i', '--input-bam', metavar='File', dest="input_bam", type=Path, required=True)
    parser.add_argument('-o', '--output-dir', metavar='Directory', dest='output_dir', type=Path,
                        required=True)
    parser.add_argument('-p', '--processes', dest='parallel_processes', type=int, default=12,
                        help='Number of logical cores that will be used to mark duplicates. Each core '
                             'processes an entire scaffold. So for the human reference genome with '
                             '24 large unique scaffolds, >23 would be ideal. >24 will only yield a minor '
                             'performance improvement. Should be at least 12 for a significant performance '
                             'gain.')
    return parser.parse_args()


def mark_duplicates_with_mate_cigar_parallel(bam_path: OneOf[str, Path], temp_dir: OneOf[str, Path],
                                             processes: int,
                                             complete_output_bam: OneOf[str, Path],  # needs to be controlled
                                             complete_output_metrics_file: OneOf[str, Path]) -> Path:
    try:
        samtools_path = Path(which('samtools'))
    except TypeError:  # TypeError: expected str, bytes or os.PathLike object, not NoneType
        raise FileNotFoundError("samtools executable not found on system path (looked through the eyes of "
                                "shutil.which()). Terminating..")  # re-raise proper exception
    try:
        picard_path = Path(which('picard'))
    except TypeError:  # TypeError: expected str, bytes or os.PathLike object, not NoneType
        raise FileNotFoundError("picard executable not found on system path (looked through the eyes of "
                                "shutil.which())")  # re-raise proper exception
    half_processes = processes // 2
    if half_processes < 1:
        half_processes = 2
    # FIX MATE INFORMATION FIRST USING ALL THE SPECIFIED CORES!
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)
    fixed_bam_path = temp_path / f"{bam_path.stem}.mateFixed.sorted.bam"
    # the chained command below uses the defined number of processes
    namesorting_tmp_dir = temp_path / f'namesorting-{bam_path.stem}'
    if namesorting_tmp_dir.is_dir():
        rmtree(namesorting_tmp_dir, ignore_errors=True)
    namesorting_tmp_dir.mkdir(parents=True, exist_ok=True)
    namesorted_temp_prefix = namesorting_tmp_dir / bam_path.stem
    # care for temporary directory of namesorting
    coordsorting_tmp_dir = temp_path / f'coordsorting-{bam_path.stem}'
    if coordsorting_tmp_dir.is_dir():
        rmtree(coordsorting_tmp_dir, ignore_errors=True)
    coordsorting_tmp_dir.mkdir(parents=True, exist_ok=True)
    coordsorted_temp_prefix = coordsorting_tmp_dir / bam_path.stem
    # assemble piped command
    mate_info_fixing_command = (f"{samtools_path} sort --threads {half_processes - 1} -n "
                                f"-T {namesorted_temp_prefix} {bam_path} | "  # SORTING: 11 threads
                                # ensure input to fixmate is name-sorted 
                                # (BAM is likely corrdinate sorted)
                                f"{samtools_path} fixmate -m - - | "  # MATE FIXING: 1 thread
                                # adds required mate tags like 'MQ' (quality) and 'MC' (CIGAR string), 
                                # if missing
                                f"{samtools_path} sort --threads {half_processes - 1} "
                                f"-l 1 "  # use lightest/fastest compression
                                f"-o {fixed_bam_path} -T {coordsorted_temp_prefix} - && "  # SORTING: 11 threads
                                f"{samtools_path} index {fixed_bam_path} {fixed_bam_path}.bai")  # BAM INDEXING: 1 thread (for portability, don't use -o option)
    # ensure input to MarkDuplicatesWithMateCigar is coordinate-sorted
    called_mate_fixing_cmd = subp.run(mate_info_fixing_command, encoding='utf-8', shell=True,
                                      stderr=subp.PIPE, stdout=subp.PIPE)
    try:
        called_mate_fixing_cmd.check_returncode()
    except subp.CalledProcessError:
        print(f"Mate information fixing child process returned with non-zero exit status.\n"
              f"This was the output: {called_mate_fixing_cmd.stdout}\n"
              f"This was the error output: {called_mate_fixing_cmd.stderr}")
        exit(3)
    # delete the temporary directories for sorting again
    if namesorting_tmp_dir.is_dir():
        rmtree(namesorting_tmp_dir, ignore_errors=True)
    if coordsorting_tmp_dir.is_dir():
        rmtree(coordsorting_tmp_dir, ignore_errors=True)

    # split reference scaffold list into approximately equal-sized lists of scaffolds (bp-wise)
    ref_seq_order_command = f"""samtools view -H {fixed_bam_path} | grep '^@SQ' | awk -F'\t' '{{
        for (i = 1; i <= NF; i++) {{ 
            if ($i ~ /^SN:/) {{ split($i, a, ":"); sn = a[2] }}
            if ($i ~ /^LN:/) {{ split($i, b, ":"); ln = b[2] }}
        }}
        print sn "," ln
    }}'"""  # names only: f"samtools view -H {fixed_bam_path} | grep '^@SQ' | cut -f2 | sed 's/^SN://'"
    called_ref_order_command = subp.run(ref_seq_order_command, shell=True, capture_output=True,
                                        encoding='utf-8')
    try:
        called_ref_order_command.check_returncode()
    except CalledProcessError:
        raise SubprocessError(f"could not determine reference sequence order. Exiting ..")
    ref_dict = {}
    reference_names = []
    _ = [(ref_dict.update({ref_name: ref_length}),
          reference_names.append(ref_name))
         for t in called_ref_order_command.stdout.strip().splitlines()
         if (ref_name := t.split(',')[0]) and (ref_length := int(t.split(',')[1]))]
    total_ref_length = sum(ref_dict.values())
    optimal_split_length = total_ref_length // processes
    # create scaffold lists for workers
    worker_scaffold_lists = [[] for _l in range(processes)]
    worker_scaffolds_cumulative_sums = [0 for _ in range(processes)]
    descending_length_tuples = sorted(list(zip(reference_names, [ref_dict[r] for r in reference_names])),
                                      key=lambda t: t[1], reverse=True)
    for scaff_idx, (scaffold_name, scaffold_length) in enumerate(descending_length_tuples):
        # check if adding to naive target is a valid choice
        if worker_scaffolds_cumulative_sums[scaff_idx % processes] + scaffold_length < optimal_split_length:  # default
            worker_scaffold_lists[scaff_idx % processes].append(scaffold_name)
            worker_scaffolds_cumulative_sums[scaff_idx % processes] += scaffold_length
        else:  # add to the list with the lowest cumulative base sum
            shortest_list_index = worker_scaffolds_cumulative_sums.index(min(worker_scaffolds_cumulative_sums))
            worker_scaffold_lists[shortest_list_index].append(scaffold_name)
            worker_scaffolds_cumulative_sums[shortest_list_index] += scaffold_length
    # create multiprocesses:
    receiving_connections, child_sending_connections = zip(*[mp.Pipe(duplex=False)   # unpack these Connections!
                                                             for _w_idx in range(parallel_processes)])
    # create worker_kwargs
    worker_kwargs = [{'bam_path': fixed_bam_path, 'temp_dir': temp_dir, 'output_bam': complete_output_bam,
                      'output_metrics_file': complete_output_metrics_file, 'picard_path': picard_path,
                      'scaffolds_to_process': worker_scaffold_lists[w_idx], 'samtools_path': samtools_path,
                      'parent_connection': child_sending_connections[w_idx]}
                     for w_idx in range(processes)]
    # create workers
    worker_processes = [mp.Process(target=mark_duplicates_with_mate_cigar, kwargs=worker_kwargs[worker_idx])
                        for worker_idx in range(processes)]
    # start all workers
    _ = deque(map(lambda p: p.start(), worker_processes), maxlen=0)
    # collect results
    received_scaffold_files = []  # (output_bam, scaffold_metrics_file)
    for r_c in receiving_connections:
        received_scaffold_files.extend(r_c.recv())  # might fail here if the storage is full:
        # srun: error: sy096: task 0: Exited with exit code 1
    # check results and terminate if one of the subprocesses failed
    if None in received_scaffold_files:
        print(f'At least one (actually {sum([rbp is None for rbp in received_scaffold_files]):,}) single '
              'scaffold duplicates marking process failed. Terminating ..')
        exit(1)
    received_scaffold_bam_paths, received_scaffold_metrics_files = zip(*received_scaffold_files)
    # samtools cat the individual scaffold BAMs in order
    # move all scaffold metrics files to a 'metrics' directory in the output_dir
    metrics_output_dir = complete_output_metrics_file.parent / f'{bam_path.stem}_scaffold_metrics'
    metrics_output_dir.mkdir(parents=True, exist_ok=True)
    for metrics_file in received_scaffold_metrics_files:  # move from scaffold temp to metrics file parent dir
        if (metrics_output_dir / metrics_file.name).is_file():  # file already exists at destination
            print(f"WARNING: mark duplicates metrics file ('{metrics_file.name}') already exists at "
                  f"destination folder '{metrics_output_dir}'. Deleting destination ..")
            (metrics_output_dir / metrics_file.name).unlink()
        move(metrics_file, metrics_output_dir)
    # TODO (restore feature broken by multiprocessing into scaffold-wise statistics): create function to
    #  aggregate the statistics in all scaffold metrics files and re-create the complete metrics file!
    # check if output BAM file exists -> delete if so!
    concatenated_bam = temp_path / complete_output_bam.name
    if concatenated_bam.is_file():
        print(f"WARNING: output BAM file already exists in temporary directory. "
              f"Deleting it before concatenating scaffold BAM files ..")
        concatenated_bam.unlink()
    # create the correct order of DupMarked scaffold BAMs
    received_order = [rcv_scaf_path.stem.split('-')[-1] for rcv_scaf_path in received_scaffold_bam_paths]
    ordered_scaffold_paths = [str(received_scaffold_bam_paths[received_order.index(order_scaff)])
                              for order_scaff in reference_names]

    concatenation_command = f"{samtools_path} cat {' '.join(ordered_scaffold_paths)} -o {concatenated_bam} "
    concatenation_subprocess = subp.run(concatenation_command, shell=True, stdout=subp.PIPE, stderr=subp.PIPE)
    try:
        concatenation_subprocess.check_returncode()
    except subp.CalledProcessError:
        print(f"concatenation of scaffold BAM files failed. Terminating ..")
        exit(1)

    # delete individual scaffold BAM files and their parent directories
    for scaffold_bam_paths in received_scaffold_bam_paths:
        rmtree(scaffold_bam_paths.parent, ignore_errors=True)  # also deletes the directory itself
    # delete the mate-fixed BAM file index
    index_files = fixed_bam_path.parent.glob(f'{fixed_bam_path.stem}*.bai')
    for idx_file in index_files:
        idx_file.unlink()
    # delete the mate-fixed BAM file and its index
    fixed_bam_path.unlink(missing_ok=True)
    # index the duplicates-marked BAM file
    index_path = Path(f"{str(concatenated_bam)}.bai")
    index_command = [f"{samtools_path}", 'index', f'{concatenated_bam}', f'{index_path}']
    indexing_subprocess = subp.run(index_command, stdout=subp.PIPE, stderr=subp.PIPE)
    try:
        indexing_subprocess.check_returncode()
    except subp.CalledProcessError:
        print(f"indexing of duplicates marked concatenated BAM file failed. Terminating ..")
        exit(1)
    # move the BAM file and the index to the target directory
    if complete_output_bam.is_file():
        print(f"WARNING: output BAM file already exists. Deleting it before "
              f"rsyncing result BAM file and index..")
        complete_output_bam.unlink()
    complete_output_bam_index_path = Path(f'{complete_output_bam}.bai')
    if complete_output_bam_index_path.is_file():
        print(f"WARNING: output BAM file index already exists. Deleting it before "
              f"rsyncing..")
        complete_output_bam_index_path.unlink()  # delete
    # synchronize to target location
    final_output_dir = complete_output_bam.parent
    final_output_dir.mkdir(exist_ok=True, parents=True)
    rsync_command = f'{rsync_path} --checksum {str(concatenated_bam)[:-1]}* {final_output_dir}/'
    # trailing slash: place inside! Use shell=True to enable shell expansion of wildcard!
    rsync_subprocess = subp.run(rsync_command, shell=True, stdout=subp.PIPE, stderr=subp.PIPE)
    try:
        rsync_subprocess.check_returncode()
    except subp.CalledProcessError:
        print(f"rsyncing duplicates marked concatenated BAM file and index failed. Terminating ..")
        exit(1)
    # delete the synchronized files (they will pop up in the )
    index_path.unlink()
    concatenated_bam.unlink()
    return temp_path  # THIS DIRECTORY MUST BE SAMPLE-SPECIFIC!


def compute_insert_size_metrics(input_bam_path: OneOf[str, Path], output_directory: OneOf[str, Path],
                                sample_id: str = None) -> subp.Popen:
    try:
        picard_path = Path(which('picard'))
    except TypeError:  # TypeError: expected str, bytes or os.PathLike object, not NoneType
        raise FileNotFoundError("picard executable not found on system path (looked through the eyes of "
                                "shutil.which())")  # re-raise proper exception
    in_bam_path = Path(input_bam_path)
    output_dir_path = Path(output_directory)  # this should be on the staging drive!
    if in_bam_path.parent == output_dir_path:
        raise ValueError(f"Output directory path must not be identical to input BAM parent directory!")
    output_dir_path.mkdir(parents=True, exist_ok=True)
    if sample_id is None:
        sample_id = in_bam_path.name.split('.markdup.bam')[0]
    metrics_output_path = output_dir_path / f'{sample_id}-insert_size_metrics.txt'
    histogram_output_path = output_dir_path / f"{sample_id}-insert_size_histogram.pdf"
    insert_sizes_tmp_dir = output_dir_path / f'{sample_id}-insert_sizes_tmp_dir'
    if insert_sizes_tmp_dir.is_dir():
        rmtree(insert_sizes_tmp_dir, ignore_errors=True)
    insert_sizes_tmp_dir.mkdir(parents=True, exist_ok=True)
    # assemble insert size metrices command:
    insert_sizes_cmd = [str(picard_path), '-Xmx15g',
                        'CollectInsertSizeMetrics',
                        f'I="{input_bam_path}"',
                        f'O="{metrics_output_path}"',
                        f'H="{histogram_output_path}"',
                        'M=0.05',
                        f'TMP_DIR={insert_sizes_tmp_dir}']
    insert_sizes_proc = subp.Popen(insert_sizes_cmd, encoding='utf-8', stderr=subp.PIPE, stdout=subp.PIPE)
    return insert_sizes_proc, insert_sizes_tmp_dir


def rsync_results_to_output_dir(source_directory: OneOf[str, Path], destination_folder: OneOf[str, Path]):
    if not destination_folder.is_dir():
        destination_folder.mkdir(parents=True, exist_ok=True)
        rmtree(destination_folder, ignore_errors=True)
    rsync_command = [str(rsync_path), '-rl', '--checksum',
                     f'{source_directory}/',  # trailing slash to synchronize all content! (insert sizes dir)
                     f'{destination_folder}/']  # trailing slash to place all content into the destination!
    rsync_subprocess = subp.run(rsync_command, stdout=subp.PIPE, stderr=subp.PIPE)
    try:
        rsync_subprocess.check_returncode()
    except subp.CalledProcessError:
        print(f"rsyncing results from duplicate marking and insert size metrics computation failed. "
              f"Terminating before deleting temporary directory '{source_directory}' ..")
        exit(1)
    rmtree(source_directory, ignore_errors=True)  # is outside the scaffold_bam_paths.parent directory!


# execution requires GATK4+ to be present with Picard MarkDuplicatesWithMateCigar
def mark_duplicates_with_mate_cigar(bam_path: OneOf[str, Path], temp_dir: OneOf[str, Path],
                                    output_bam: OneOf[str, Path], output_metrics_file: OneOf[str, Path],
                                    scaffolds_to_process: List[str], parent_connection: mp.Pipe,
                                    picard_path: OneOf[str, Path], samtools_path: OneOf[str, Path]):
    created_files = []
    # create duplicates-marked scaffold BAMs
    for scaffold in scaffolds_to_process:
        scaffold_temp_dir = temp_dir / f'{scaffold}-{randint(0, 9999999):07}/tmp'
        scaffold_temp_dir.mkdir(parents=True, exist_ok=True)
        scaffold_output_bam = scaffold_temp_dir.parent / f'{output_bam.stem}-{scaffold}.bam'
        scaffold_metrics_file = scaffold_temp_dir.parent / (f"{output_metrics_file.stem}-{scaffold}."
                                                            f"{output_metrics_file.name.split('.')[-1]}")
        # samtools filters out unaligned reads before marking duplicates
        # (command terminates otherwise if unaligned reads are present)
        # -F 12 means remove any read pair for which either the read or the mate is unmapped
        # WARNING - this command requires a correct software environment to be active!
        # (e.g., TSO500ichorCNA of Benjamin on MedBioNode)

        # the chained command below uses 2 processes -> scaffold selection & duplicates marking
        paired_end_dup_marking_cmd = (f"{samtools_path} view -h -F 12 --uncompressed {bam_path} {scaffold} | "
                                      # 24 of these processes will use max of 24x2 = 48 Gb
                                      f"{picard_path} -Xmx2g MarkDuplicatesWithMateCigar "  # 1 thread
                                      f"INPUT=/dev/stdin "
                                      f"TMP_DIR={scaffold_temp_dir}/tmp "
                                      f"OUTPUT={scaffold_output_bam} "
                                      f"METRICS_FILE={scaffold_metrics_file} "
                                      f"ASSUME_SORTED=true "
                                      # why can I use /dev/stdin and why does it not interfere between 
                                      # different processes?
                                      # ++++++++++++++++++++
                                      # /dev/stdin is a symbolic link to /proc/self/fd/0, and THUS process specific!
                                      # /proc/self is a symbolic link only seen by your running process to its 
                                      # process-id; The /proc filesystem is a virtual (not real) filesystem which can 
                                      # show a different view to each process.
                                      # --MINIMUM_DISTANCE <Integer>
                                      #   The minimum distance to buffer records to account for clipping on 
                                      #   the 5' end of the records. For a given alignment, this parameter 
                                      #   controls the width of the window to search for duplicates of that 
                                      #   alignment. Due to 5' read clipping, duplicates do not necessarily
                                      #   have the same 5' alignment coordinates, so the algorithm needs to 
                                      #   search around the neighborhood. For single end sequencing data, the 
                                      #   neighborhood is only determined by the amount of clipping (assuming 
                                      #   no split reads), thus setting MINIMUM_DISTANCE to twice the 
                                      #   sequencing read length should be sufficient. For paired end 
                                      #   sequencing, the neighborhood is also determined by the fragment 
                                      #   insert size, so you may want to set MINIMUM_DISTANCE to something 
                                      #   like twice the 99.5% percentile of the fragment insert size 
                                      #   distribution (see CollectInsertSizeMetrics). Or you can set this 
                                      #   number to -1 to use either a) twice the first read's read length, 
                                      #   or b) 100, whichever is smaller. Note that the larger the window, 
                                      #   the greater the RAM requirements, so you could run into performance 
                                      #   limitations if you use a value that is unnecessarily large.
                                      #   Default value: -1.
                                      #   99th percentile insert size for samples of run '250627_TSO500_Onco':
                                      #   395 bp, 409 bp, 427 bp, 437 bp, 475 bp, 607 bp, 631 bp
                                      "MINIMUM_DISTANCE=750 "
                                      "COMPRESSION_LEVEL=2 "  # default is 2
                                      "CREATE_INDEX=true "  # Creates '{file_path.stem}.bai' index file
                                      # --MAX_RECORDS_IN_RAM <Integer>
                                      #   When writing files that need to be sorted, this will specify the 
                                      #   number of records stored in RAM before spilling to disk. Increasing 
                                      #   this number reduces the number of file handles needed to sort the 
                                      #   file, and increases the amount of RAM needed.
                                      #   Default value: 500000.
                                      "MAX_RECORDS_IN_RAM=1000000 "
                                      # Illumina patterned flowcell:
                                      # --OPTICAL_DUPLICATE_PIXEL_DISTANCE <Integer>
                                      #   The maximum offset between two duplicate clusters to consider them 
                                      #   optical duplicates. The default is appropriate for unpatterned 
                                      #   versions of the Illumina platform. For the patterned flowcell 
                                      #   models, 2500 is more appropriate. For other platforms and models, 
                                      #   users should experiment to find what works best.
                                      #   Default value: 100.
                                      "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 "
                                      "REMOVE_DUPLICATES=false")  # make it explicit -> we only tag
        # --REMOVE_DUPLICATES <Boolean> If true do not write duplicates to the output file instead of
        #   writing them with appropriate flags set.
        #   Default value: false.
        #   Possible values: {true, false}
        #
        # detected about 15% duplication after all in sample ~10 GB input BAM file:
        # (BAM file shrinks to 8.4 GB after fixing mate information)
        # '/home/isilon/HumGenTempData/for_Raul_TSO500/Data_for_ichor/250611_TSO500_Onco-RECREATIONNEWSCRIPT/
        #  L186_1-ONC_tumor.mateFixed.sorted.bam'
        paired_end_dup_marking_cmd = subp.run(paired_end_dup_marking_cmd, encoding='utf-8', shell=True,
                                              stderr=subp.PIPE, stdout=subp.PIPE)
        try:
            paired_end_dup_marking_cmd.check_returncode()
            created_files.append((scaffold_output_bam, scaffold_metrics_file))
        except subp.CalledProcessError:
            print(f"Duplicates marking child process returned with non-zero exit status for scaffold "
                  f"'{scaffold}'.\n"
                  f"This was the output: {paired_end_dup_marking_cmd.stdout}\n"
                  f"This was the error output: {paired_end_dup_marking_cmd.stderr}")
            created_files.append(None)
            break  # one failure is enough to quit
    if None in created_files:
        parent_connection.send(None)
    else:
        parent_connection.send(created_files)  # might send an empty list if no scaffolds were received
    parent_connection.close()


def wait_for_insert_metrics(process_handle: subp.Popen, metrics_output_dir: OneOf[str, Path],
                            delete_temp_dir: OneOf[str, Path], sample_id: str, max_wait_minutes: int = 180):
    wait_time = 0
    waited_hours = waited_minutes = waited_seconds = 0
    while process_handle.poll() is None:
        print("waiting for insert metrics to finish..")
        sleep(10)
        wait_time += 10
        waited_hours, remainder = divmod(wait_time, 3600)
        waited_minutes, waited_seconds = divmod(remainder, 60)
        if wait_time // 60 >= max_wait_minutes:
            waited_hours, remainder = divmod(wait_time, 3600)
            waited_minutes, waited_seconds = divmod(remainder, 60)
            process_handle.terminate()
            process_handle.wait() # close the process!
            if process_handle.stdout:
                process_handle.stdout.close()
            if process_handle.stderr:
                process_handle.stderr.close()
            raise TimeoutError(f"waited now {waited_hours:02}:{waited_minutes:02}:{waited_seconds:02} "
                               "(hh:mm:ss) for insert size metrics to finish. I don't wanna wait anymore!")
    # read stdout und stderr:
    process_stdout, process_stderr = process_handle.communicate()
    # process finished; move on
    if wait_time:
        print(f"INFO - waited {waited_hours:02}:{waited_minutes:02}:{waited_seconds:02} (hh:mm:ss) "
              f"for insert size metrics to finish.")
    if process_handle.returncode:
        print(f"Insert sizes computation child process returned with non-zero exit.\n"
              f"This was the output: {process_stdout}\n"
              f"This was the error output: {process_stderr}")
        exit(1)
    output_dir_path = Path(metrics_output_dir)
    # create a subdirectory containing all the insert size metrics
    insert_size_files_dir_path = output_dir_path / 'insert_sizes'
    insert_size_files_dir_path.mkdir(parents=True, exist_ok=True)  # may have been created by another process
    content_pattern = f'{sample_id}-insert_size_*'
    for file_objct in output_dir_path.glob(content_pattern):
        if file_objct.is_file():
            move(file_objct, f'{insert_size_files_dir_path}/{Path(file_objct).name}')
        elif file_objct.is_dir():
            move(file_objct, f'{insert_size_files_dir_path}/{Path(file_objct).name}/')
    # delete obsolete insert size metrics tmp directory
    if delete_temp_dir.is_dir():
        rmtree(delete_temp_dir)


if __name__ == '__main__':
    # get commandline arguments
    in_args = get_cmdline_args()
    input_bam = Path(str(in_args.input_bam).strip('"').strip("'"))
    sample_name = Path(input_bam).name.split('.markdup.bam')[0]
    output_dir = in_args.output_dir
    temp_dir = in_args.temp_dir
    parallel_processes = in_args.parallel_processes
    if parallel_processes > 24:
        print(f"WARNING - {parallel_processes:,} parallel processes were defined but maximum allowed is 32. "
              f"Limiting to 32.")
        parallel_processes = 24
    elif parallel_processes < 2:
        print(f"WARNING - {parallel_processes:,} parallel processes were defined but at least 2 should be "
              f"used. Setting to 2.")
        parallel_processes = 2
    # construct missing parameters
    metrics_file = Path(output_dir) / f"{input_bam.stem}.markdup.metrics"
    output_bam_path = Path(output_dir) / f"{input_bam.stem}.markdup.bam"
    output_dir.mkdir(parents=True, exist_ok=True)
    # STEP1a: asynchronously compute insert size metrics from duplicates-marked BAM file
    insert_metrics_process, insert_szs_tmp_path = compute_insert_size_metrics(
        input_bam_path=input_bam, sample_id=sample_name, output_directory=temp_dir)
    # STEP1b: function below runs samtools fixmate and then executes picard MarkDuplicatesWithMateCigar
    rsync_this_dirs_contents = mark_duplicates_with_mate_cigar_parallel(
        bam_path=input_bam, temp_dir=temp_dir, processes=parallel_processes,
        complete_output_bam=output_bam_path, complete_output_metrics_file=metrics_file)
    # wait for asynchronous insert sizes process to finish
    wait_for_insert_metrics(
        process_handle=insert_metrics_process, sample_id=sample_name,
        metrics_output_dir=temp_dir, delete_temp_dir=insert_szs_tmp_path)  # MUST be 'output_directory' from compute_insert_size_metrics() !
    # move results from temp to output directory
    rsync_results_to_output_dir(source_directory=rsync_this_dirs_contents,
                                destination_folder=output_dir)
