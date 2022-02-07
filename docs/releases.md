
# Release binaries

Binary executables for releases of the Rust client are available for Linux and MacOS.

<div class="termy">

Configure environmental variables.

```console
VERSION=0.5.0
GITHUB=https://github.com/esteinig/sketchy/releases/download
```

Download release binaries (Linux).

```console
$ TAR=sketchy-${VERSION}-x86_64-unknown-linux-musl.tar.gz
$ wget ${GITHUB}/${VERSION}/${TAR}
---> 100%
$ tar xf $TAR
---> 100%
```

Download release binaries (MacOS).

```console
$ TAR=sketchy-${VERSION}-x86_64-apple-darwin.tar.gz
$ wget ${GITHUB}/${VERSION}/${TAR}
---> 100%
$ tar xf $TAR
---> 100%
```

</div>