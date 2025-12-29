FROM continuumio/miniconda3:latest

WORKDIR /mnt/

COPY ichorCNA.yaml /tmp/ichorCNA.yaml

RUN conda env create -f /tmp/ichorCNA.yaml && conda clean -afy

SHELL ["/bin/bash", "-lc"]

RUN echo "conda activate TSO500ichorCNA" >> ~/.bashrc

ENV ICHOR_SCRIPT="tso500_dragen_pipeline/sandbox/ichorCNA_docker/tso500offtarget_ichorcna/drv_TSO500_offtarget_ichorCNA_parallel.sh" \
    DATA_DIR="tso500_dragen_pipeline/sandbox/ichorCNA_docker/data" \
    ICHOR_DIR="tso500_dragen_pipeline/ichorCNA"
