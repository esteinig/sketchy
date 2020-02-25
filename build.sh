mkdir -p $PREFIX/bin

cargo build --release
mv target/release/sketchy $PREFIX/bin/sketchy-rs

$PYTHON -m pip install .
