mkdir -p $PREFIX/bin

cargo build --release
mv target/release/sketchy-rs $PREFIX/bin/sketchy-rs

$PYTHON -m pip install .
