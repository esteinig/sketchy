Bootstrap: docker
From: continuumio/miniconda3

%labels
    name sketchy
    version 0.4.4
    author esteinig

%post
    export PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
    export SKETCHY_PATH=/sketchy

    apt-get update && apt-get install curl build-essential -y

    mkdir /rust
    export CARGO_HOME=/rust/.cargo
    export RUSTUP_HOME=/rust/.rustup

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

    cd /sketchy_build \
        && cargo build --release \
        && ls /sketchy_build/target/release \
        && mv /sketchy_build/target/release/sketchy-rs /bin/sketchy-rs

    /opt/conda/bin/conda install -c conda-forge -c bioconda -c esteinig --yes \
        pysam mash=2.2 psutil nanoq \
        && /opt/conda/bin/conda clean -a \
        && find /opt/conda/ -follow -type f -name '*.a' -delete \
        && find /opt/conda/ -follow -type f -name '*.pyc' -delete

    pip install /sketchy_build
    sketchy pull --full --path $SKETCHY_PATH

%environment
    export PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
    export SKETCHY_PATH=/sketchy

%files
    . /sketchy_build

%test
    export PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
    sketchy --help
    mash -h
    sketchy-rs -h

