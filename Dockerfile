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

RUN conda install -c conda-forge -c bioconda -c esteinig --yes pysam mash=$MASH_VERSION

RUN mkdir /rust && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Latest Rust client and default sketches
RUN cargo install sketchy-rs
RUN mkdir /sketches && sketchy get --outdir /sketches

# Latest Python utility client
RUN pip install git+https://github.com/esteinig/sketchy

