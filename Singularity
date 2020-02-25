Bootstrap: docker
From: continuumio/miniconda3

%labels
    name sketchy
    version 0.4.2
    author esteinig

%post
    export PATH=/opt/conda/bin:/rust/.cargo/bin:$PATH
    export SKETCHY_PATH=/sketchy

    apt-get update && apt-get install curl build-essential -y

    mkdir /rust
    export CARGO_HOME=/rust/.cargo
    export RUSTUP_HOME=/rust/.rustup

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

    /opt/conda/bin/conda install -c conda-forge -c bioconda -c esteinig --yes \
        pysam mash=2.2 psutil pf-core \
        && /opt/conda/bin/conda clean -a \
        && find /opt/conda/ -follow -type f -name '*.a' -delete \
        && find /opt/conda/ -follow -type f -name '*.pyc' -delete

    cd /sketchy_build \
        && cargo build --release \
        && mv /sketchy_build/target/release/sketchy /bin/sketchy-rs &&
        pip install .

    sketchy pull --path $SKETCHY_PATH

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

