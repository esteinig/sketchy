FROM continuumio/miniconda3

LABEL name="sketchy"
LABEL version="0.5.0"
LABEL author="esteinig"

RUN apt-get update && apt-get install curl build-essential -y

ENV SKETCHY_VERSION=0.5.0
ENV MASH_VERSION=2.2.2
ENV CONDA_DIR=/opt/conda
ENV SKETCHY_PATH=/sketchy

ENV PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
ENV CARGO_HOME=/rust/.cargo
ENV RUSTUP_HOME=/rust/.rustup

RUN conda install -c conda-forge -c bioconda -c esteinig --yes \
    pysam mash=$MASH_VERSION psutil nanoq

RUN mkdir /rust && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN cargo install sketchy-rs

RUN wget https://github.com/esteinig/sketchy/archive/v${SKETCHY_VERSION}.tar.gz
RUN tar -xzf v${SKETCHY_VERSION}.tar.gz && rm v${SKETCHY_VERSION}.tar.gz
RUN pip install /sketchy-${SKETCHY_VERSION}



