FROM continuumio/miniconda3

LABEL name="sketchy"
LABEL version="0.4.0"
LABEL author="esteinig"

RUN apt-get update && apt-get install curl build-essential -y

ENV MASH_VERSION="2.2"
ENV CONDA_DIR="/opt/conda"
ENV SKETCHY_PATH="/sketchy"
ENV PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
ENV CARGO_HOME=/rust/.cargo
ENV RUSTUP_HOME=/rust/.rustup

RUN conda install -c conda-forge -c bioconda -c esteinig --yes \
    pysam mash=$MASH_VERSION psutil pf-core \
    && conda clean -a \
    && find $CONDA_DIR -follow -type f -name '*.a' -delete \
    && find $CONDA_DIR -follow -type f -name '*.pyc' -delete

RUN mkdir /rust && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

ADD . /sketchy_build
WORKDIR /sketchy_build

# Sketchy Rust library
RUN cargo build --release && mv /sketchy_build/target/release/sketchy /bin/sketchy-rs
# Sketchy Python package and data working directory
RUN pip install /sketchy_build && mkdir /data
# Sketchy database pull
RUN sketchy pull --path $SKETCHY_PATH

# /data workdir for easy bindmounts
WORKDIR /data
ENTRYPOINT ["sketchy"]
CMD []


