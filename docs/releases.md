
# Release binaries

Binary executables for releases of the Rust client are available for Linux and MacOS.

Configure environmental variables.

```
VERSION=0.6.0
GITHUB=https://github.com/esteinig/sketchy/releases/download
```

Download release binaries (Linux).

```
TAR=sketchy-${VERSION}-x86_64-unknown-linux-musl.tar.gz
wget ${GITHUB}/${VERSION}/${TAR}
100%
tar xf $TAR
```

Download release binaries (MacOS).

```
TAR=sketchy-${VERSION}-x86_64-apple-darwin.tar.gz
wget ${GITHUB}/${VERSION}/${TAR}
tar xf $TAR
```

