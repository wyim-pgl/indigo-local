# indigo-local

Local command-line setup of [Indigo](https://github.com/gear-genomics/indigo) (gear-genomics) for detecting CRISPR/TALEN-induced SNVs and InDels in Sanger chromatograms, with two additions:

- **`indigo_report.py`** — renders the tracy JSON output as a standalone HTML report that looks like the Indigo web app (variants table, trace viewer, decomposition plot, alignments).
- **`indigo_primer.py`** — primer-based trimming: finds the PCR primers in the basecalls and sets the left/right trim sizes before running Indigo.

The original Indigo README is kept in [README.original.md](README.original.md). Indigo and Tracy are developed by Tobias Rausch (EMBL) and licensed under GPL-3.0 / BSD-3 respectively.

## Requirements

- [Tracy](https://github.com/gear-genomics/tracy) (`tracy` in `PATH`)
- R with `gplots`, `ggplot2`, `reshape2`
- Python 3 (standard library only)

## Installation

### Option A: bioconda (easiest, not tested in this setup)

```bash
micromamba create -n indigo -c conda-forge -c bioconda tracy r-base r-gplots r-ggplot2 r-reshape2 python
micromamba activate indigo
git clone https://github.com/wyim-pgl/indigo-local.git
```

### Option B: build Tracy from source inside a conda environment (tested on WSL2 Ubuntu)

```bash
# system build tools
sudo apt update
sudo apt install -y build-essential autoconf automake cmake \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libboost-all-dev

# conda dependencies
micromamba activate <env>
micromamba install -c conda-forge -c bioconda htslib boost-cpp cmake autoconf automake make zlib bzip2 xz
micromamba install -c conda-forge r-gplots r-ggplot2 r-reshape2

# build tracy against the conda htslib/boost
git clone --recursive https://github.com/gear-genomics/tracy.git
cd tracy
export CMAKE_POLICY_VERSION_MINIMUM=3.5
make HTSLIBINCDIR=$CONDA_PREFIX/include HTSLIBLIBDIR=$CONDA_PREFIX/lib \
     BOOSTINCDIR=$CONDA_PREFIX/include BOOSTLIBDIR=$CONDA_PREFIX/lib all
./src/tracy -h
export PATH=$PWD/src:$PATH
```

#### Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `autoreconf: not found` | build tools missing | install the apt packages above |
| `C compiler cannot create executables` (htslib configure) | an active conda env sets `CC`/`LDFLAGS`; tracy's `-lhts -lboost_*` leak into htslib's configure | build against the conda htslib/boost with `HTSLIBINCDIR`/`HTSLIBLIBDIR`/`BOOSTINCDIR`/`BOOSTLIBDIR` |
| `Compatibility with CMake < 3.5 has been removed` | CMake 4.x vs. bundled sdsl | `export CMAKE_POLICY_VERSION_MINIMUM=3.5` |
| `cannot find -lz` | zlib missing in the conda env | `micromamba install -c conda-forge zlib bzip2 xz` |
| `DYLD_LIBRARY_PATH` has no effect | macOS-only variable | not needed on Linux; the binary has an rpath |

Do not run two builds in the same directory at the same time (`clean.sh` deletes the other build's files).

## Usage

### Indigo (fixed trim sizes)

```bash
./indigo.sh <experiment.ab1> <reference.fa|.fa.gz|wildtype.ab1> <ltrim> <rtrim> <pratio> <outprefix>
./indigo.sh server/sample.abi server/sample.fa 50 50 0.33 out
```

| Argument | tracy decompose option | Web default |
|---|---|---|
| `ltrim` | `-q` / `--trimLeft` | 50 |
| `rtrim` | `-u` / `--trimRight` | 50 |
| `pratio` | `-p` / `--pratio` | 0.33 (33 on the web slider) |

Outputs: `<outprefix>.pdf`, `.html` (web-like report), `.decomp`, `.bcf`, `.align1-3`, `.json`.

### Indigo with primer trimming

```bash
python3 indigo_primer.py <experiment.ab1> <reference.fa> -f <fwd primer> -r <rev primer> -o <outprefix>
python3 indigo_primer.py sample.ab1 ref.fa -f ACGT... -r TGCA... --dry-run   # only report trim sizes
```

| Option | Meaning | Default |
|---|---|---|
| `-f`, `-r` | primer sequences 5'→3' (one is enough) | — |
| `-m` | max. mismatches (IUPAC codes allowed) | 2 |
| `--pad` | extra bases trimmed beyond a primer | 0 |
| `--ltrim`, `--rtrim` | trim size used when no primer is found on that side | 50, 50 |
| `-p` | peak ratio | 0.33 |

Primers and their reverse complements are searched in the primary basecalls. A hit in the 5' half sets the left trim to the primer end, a hit in the 3' half sets the right trim from the primer start. Note that in a typical Sanger read the sequencing primer itself is not present, so usually only the opposite PCR primer is found (near the 3' end); the other side falls back to `--ltrim`/`--rtrim`.

### HTML report only

```bash
python3 indigo_report.py <outprefix>.json [report.html]
```

The report loads Plotly and Bootstrap from a CDN, so an internet connection is needed to view it.

## Example

[`example/`](example/) contains the result for the bundled sample (`server/sample.abi` vs. `server/sample.fa`, trims 50/50, pratio 0.33). The decomposition has its minimum at **−7 bp** — a heterozygous 7 bp deletion (`chrA:1226 GGAGGGTC>G`) with allele fractions of ~0.49 each.

## License

GPL-3.0, as the original Indigo (see [LICENSE](LICENSE)).
