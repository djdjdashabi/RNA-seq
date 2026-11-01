# RNA editing workflow

Run Section 7 from `Supplemental_files/scripts/`. The workflow converts the six native caller outputs to a shared site format, annotates candidates, compares sets, and produces the manuscript plots.

## Files and paths

Keep `scripts/`, `annotation/`, `results/`, `reference/`, `bam/`, and `logs/` inside `Supplemental_files/`. Configuration paths are relative to `scripts/`. No sibling `essay` directory is needed.

Place these inputs in `annotation/`:

- `Homo_sapiens.GRCh38.110.gtf`
- `GRCh38_RepeatMasker_Alu.sorted.bed`
- `Homo_sapiens.gene_info.gz`
- `gene2go.gz`

The six caller files belong in `results/`: `redi_69.out`, `redi_70.out`, `sprint_69.out`, `sprint_70.out`, `jacusa_69.out`, and `jacusa_70.out`.

The NCBI sources downloaded on 31 October 2026 are identified below. Preserve these files when reproducing this analysis; newer downloads can change GO results.

| Source | SHA-256 |
| --- | --- |
| Homo_sapiens.gene_info.gz | d4c215ec6ae6c5922611ce3cfea9b5f899b547ee96d820b0a71b769d1c345b95 |
| gene2go.gz | 9eb2c43ca7c5b938560ae03301d676871e65909855ffcba8ddad67c9568a75e1 |

`build_gene2go.py` retains tax_id 9606 records and writes the human tables and `gene2go_build_manifest.tsv` in `annotation/`. Snakemake manages this step; a separate manual build is unnecessary.

## Run the analysis

From the project root:

```bash
conda activate snakemake_env
cd scripts
snakemake -n -j 8
snakemake -j 8
```

The first command previews the jobs. The second runs missing or outdated steps. To rebuild all downstream outputs from the existing caller files, use `snakemake --forceall -j 8`. Run `snakemake -n -j 8` afterward to confirm that targets are current.

The workflow overwrites its derived outputs. Keep the native caller files unchanged. Standardization initializes `strand_consistent`; gene annotation fills it from gene-strand compatibility. The displayed GO analysis uses `strand_mode=all`.

Script `figure1` through `figure4` correspond to manuscript Figures 2 through 5. Figure 1 is supplied separately. GO bubble selection includes terms significant in at least one sample-caller analysis; other displayed combinations need not be significant.

## Environments

The definitions below were exported from the installed Linux environments on 16 October 2026. They retain package versions and builds, without machine-specific prefixes. Fresh installation has not been tested; pip entries have version pins but no artifact hashes.

Use Linux x86-64 or Ubuntu WSL2 with Conda initialized. Create only missing environments. Keep SPRINT in its Python 2.7 environment with SAMtools 1.6. Stage SPRINT 0.1.8 in `SPRINT/` and `JACUSA_v2.1.16.jar` in `JACUSA2/` for the caller steps. These third-party tools are separate from the downstream scripts.

## Exported environment definitions


### rnaseq

Save the YAML block as environment.rnaseq.yml, then run:

```bash
conda env create -f environment.rnaseq.yml
```

```yaml
name: rnaseq
channels:
  - conda-forge
  - bioconda
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - _r-mutex=1.0.1=anacondar_1
  - alsa-lib=1.2.11=hd590300_1
  - binutils_impl_linux-64=2.46.1=bootstrap_h59bd682_2
  - bwidget=1.10.1=hfc2019e_3
  - bzip2=1.0.8=hda65f42_10
  - c-ares=1.34.8=hebe6cf0_2
  - ca-certificates=2026.7.22=hbd8a1cb_0
  - cairo=1.18.0=h3faef2a_0
  - curl=8.8.0=he654da7_1
  - expat=2.8.1=hecca717_1
  - fastp=0.23.4=h125f33a_5
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.14.2=h14ed4e7_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - freetype=2.12.1=h267a509_2
  - fribidi=1.0.16=hb03c661_1
  - gatk4=4.6.2.0=py310hdfd78af_1
  - gcc_impl_linux-64=16.2.0=h176d5d0_4
  - gfortran_impl_linux-64=16.2.0=hf930d89_4
  - giflib=5.2.2=ha257d8a_1
  - graphite2=1.3.15=h54a6638_1
  - gxx_impl_linux-64=16.2.0=h0d273dc_4
  - harfbuzz=9.0.0=hfac3d4d_0
  - hisat2=2.2.1=h503566f_8
  - htslib=1.21=h5efdd21_0
  - icu=73.2=h59595ed_0
  - isa-l=2.32.1=hb03c661_0
  - kernel-headers_linux-64=5.14.0=he073ed8_3
  - keyutils=1.6.3=h7cc23a3_1
  - krb5=1.21.3=h659f571_0
  - lcms2=2.16=hb7c19ff_0
  - ld_impl_linux-64=2.46.1=bootstrap_ha15bf96_2
  - lerc=4.2.0=hdb68285_0
  - libblas=3.11.0=10_h4a7cf45_openblas
  - libcups=2.3.3=h4637d8d_4
  - libcurl=8.8.0=hca28451_1
  - libdeflate=1.25=hd45a770_1
  - libedit=3.1.20250104=pl5321h373387f_1
  - libev=4.33=h280c20c_3
  - libexpat=2.8.1=hecca717_1
  - libffi=3.7.0=h81df57d_1
  - libgcc=16.2.0=ha9f2e26_4
  - libgcc-devel_linux-64=16.2.0=he3ce08f_104
  - libgcc-ng=16.2.0=h69a702a_4
  - libgfortran=16.2.0=h69a702a_4
  - libgfortran-ng=16.2.0=h69a702a_4
  - libgfortran5=16.2.0=h6b99dfc_4
  - libglib=2.80.2=hf974151_0
  - libgomp=16.2.0=he0feb66_4
  - libiconv=1.18=h0cb94f2_3
  - libjpeg-turbo=3.2.0=hb03c661_1
  - liblapack=3.11.0=10_h47877c9_openblas
  - liblzma=5.8.3=hb03c661_1
  - liblzma-devel=5.8.3=hb03c661_1
  - libnghttp2=1.58.0=h47da74e_1
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.34=pthreads_hcf972fe_1
  - libpng=1.6.43=h2797004_0
  - libsanitizer=16.2.0=h3048135_4
  - libsqlite=3.46.0=hde9e2c9_0
  - libssh2=1.11.0=h0841786_0
  - libstdcxx=16.2.0=h934c35e_4
  - libstdcxx-devel_linux-64=16.2.0=h86e191b_104
  - libstdcxx-ng=16.2.0=hdf11a46_4
  - libtiff=4.6.0=h1dd3fc0_3
  - libuuid=2.42.2=h5347b49_0
  - libwebp-base=1.6.0=hd42ef1d_1
  - libxcb=1.15=h0b41bf4_0
  - libxcrypt=4.4.38=h280c20c_0
  - libxml2=2.12.7=hc051c1a_1
  - libzlib=1.2.13=h4ab18f5_6
  - make=4.4.1=hb03c661_3
  - ncbi-vdb=3.4.1=hd63eeec_0
  - ncurses=6.6=hdb14827_1
  - openjdk=22.0.1=hb622114_0
  - openssl=3.6.4=h781a0a9_0
  - ossuuid=1.6.2=h5888daf_1001
  - packaging=26.3=pyhc364b38_0
  - pango=1.54.0=h4c5309f_1
  - pcre2=10.43=hcad00b1_0
  - perl=5.22.2.1=0
  - perl-threaded=5.32.1=hdfd78af_1
  - perl-uri=1.71=0
  - perl-xml-libxml=2.0124=0
  - perl-xml-namespacesupport=1.11=0
  - perl-xml-sax=0.99=0
  - perl-xml-sax-base=1.08=0
  - picard=3.0.0=hdfd78af_1
  - pip=26.2.1=pyh8b19718_0
  - pixman=0.46.4=h54a6638_3
  - pthread-stubs=0.4=hb03c661_1003
  - python=3.10.14=hd12c33a_0_cpython
  - r-base=4.4.0=h019f4a6_1
  - readline=8.3=hd6e31c0_1
  - samtools=1.17=hd87286a_2
  - sed=4.10=h19d0853_0
  - setuptools=84.0.0=pyh332efcf_0
  - sra-tools=3.2.0=h4304569_0
  - sysroot_linux-64=2.34=h087de78_3
  - tk=8.6.13=noxft_h4845f30_101
  - tktable=2.10=h8d826fa_7
  - tzdata=2026c=h151e31d_0
  - wheel=0.48.0=pyhd8ed1ab_0
  - xorg-fixesproto=5.0=hb9d3cd8_1003
  - xorg-inputproto=2.3.2=hb9d3cd8_1003
  - xorg-kbproto=1.0.7=hb9d3cd8_1003
  - xorg-libice=1.1.2=h280c20c_0
  - xorg-libsm=1.2.6=h0d788c3_1
  - xorg-libx11=1.8.9=h8ee46fc_0
  - xorg-libxau=1.0.12=hb03c661_2
  - xorg-libxdmcp=1.1.5=hb03c661_2
  - xorg-libxext=1.3.4=h0b41bf4_2
  - xorg-libxfixes=5.0.3=h7f98852_1004
  - xorg-libxi=1.7.10=h4bc722e_1
  - xorg-libxrender=0.9.11=hd590300_0
  - xorg-libxt=1.3.0=hd590300_1
  - xorg-libxtst=1.2.5=h4bc722e_0
  - xorg-recordproto=1.14.2=hb9d3cd8_1003
  - xorg-renderproto=0.11.1=hb9d3cd8_1003
  - xorg-xextproto=7.3.0=hb9d3cd8_1004
  - xorg-xproto=7.0.31=hb9d3cd8_1008
  - xz=5.8.3=ha02ee65_1
  - xz-gpl-tools=5.8.3=ha02ee65_1
  - xz-tools=5.8.3=hb03c661_1
  - zlib=1.2.13=h4ab18f5_6
  - zstd=1.5.6=ha6fb4c9_0
  - pip:
      - pysam==0.24.0
      - reditools3==3.7
```

### sprint_py2

Save the YAML block as environment.sprint_py2.yml, then run:

```bash
conda env create -f environment.sprint_py2.yml
```

```yaml
name: sprint_py2
channels:
  - conda-forge
  - bioconda
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - bzip2=1.0.8=hda65f42_10
  - c-ares=1.34.8=hebe6cf0_2
  - ca-certificates=2026.7.22=hbd8a1cb_0
  - certifi=2019.11.28=py27h8c360ce_1
  - curl=7.87.0=h6312ad2_0
  - keyutils=1.6.3=h7cc23a3_1
  - krb5=1.20.1=hf9c8cef_0
  - ld_impl_linux-64=2.46.1=bootstrap_ha15bf96_2
  - libcurl=7.87.0=h6312ad2_0
  - libedit=3.1.20250104=pl5321h373387f_1
  - libev=4.33=h280c20c_3
  - libffi=3.2.1=he1b5a44_1007
  - libgcc=16.2.0=ha9f2e26_4
  - libgcc-ng=16.2.0=h69a702a_4
  - libgomp=16.2.0=he0feb66_4
  - libnghttp2=1.51.0=hdcd2b5c_0
  - libsqlite=3.46.0=hde9e2c9_0
  - libssh2=1.10.0=haa6b8db_3
  - libstdcxx=16.2.0=h934c35e_4
  - libstdcxx-ng=16.2.0=hdf11a46_4
  - libzlib=1.2.13=h4ab18f5_6
  - ncurses=6.6=hdb14827_1
  - openssl=1.1.1w=hd590300_0
  - pip=20.1.1=pyh9f0ad1d_0
  - python=2.7.15=h5a48372_1011_cpython
  - python_abi=2.7=1_cp27mu
  - readline=8.3=hd6e31c0_1
  - samtools=1.6=hcd7b337_9
  - setuptools=44.0.0=py27_0
  - sqlite=3.46.0=h6d4b2fc_0
  - tk=8.6.13=noxft_h4845f30_101
  - wheel=0.37.1=pyhd8ed1ab_0
  - xz=5.2.6=h166bdaf_0
  - zlib=1.2.13=h4ab18f5_6
```

### jacusa

Save the YAML block as environment.jacusa.yml, then run:

```bash
conda env create -f environment.jacusa.yml
```

```yaml
name: jacusa
channels:
  - conda-forge
  - bioconda
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - alsa-lib=1.2.16.1=h7cc23a3_1
  - bzip2=1.0.8=hda65f42_10
  - c-ares=1.34.8=hebe6cf0_2
  - ca-certificates=2026.7.22=hbd8a1cb_0
  - cairo=1.18.4=h3c89d7e_3
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.18.3=h4db4eae_1
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - giflib=6.1.3=h280c20c_1
  - graphite2=1.3.15=h54a6638_1
  - htslib=1.24=ha79157c_0
  - icu=78.3=py310h44b86e0_2
  - keyutils=1.6.3=h7cc23a3_1
  - krb5=1.22.2=hbc21106_2
  - lcms2=2.19.1=h9073bf1_2
  - lerc=4.2.0=hdb68285_0
  - libcups=2.3.3=h7a8fb5f_6
  - libcurl=8.21.0=ha042cf0_5
  - libdeflate=1.25=hd45a770_1
  - libedit=3.1.20250104=pl5321h373387f_1
  - libev=4.33=h280c20c_3
  - libexpat=2.8.1=hecca717_1
  - libffi=3.7.0=h81df57d_1
  - libfreetype=2.14.3=ha770c72_2
  - libfreetype6=2.14.3=h5e6c136_2
  - libgcc=16.2.0=ha9f2e26_4
  - libglib=2.88.3=he503a2a_2
  - libgomp=16.2.0=he0feb66_4
  - libharfbuzz=14.4.0=h23af247_0
  - libiconv=1.18=h0cb94f2_3
  - libjpeg-turbo=3.2.0=hb03c661_1
  - liblzma=5.8.3=hb03c661_1
  - libnghttp2=1.68.1=h74cf4be_1
  - libpng=1.6.58=h922cc85_1
  - libpsl=0.23.1=hd9e3e90_1
  - libssh2=1.11.1=h6154650_1
  - libstdcxx=16.2.0=h934c35e_4
  - libtiff=4.7.2=hcc2c06a_1
  - libuuid=2.42.2=h5347b49_0
  - libwebp-base=1.6.0=hd42ef1d_1
  - libxcb=1.17.0=hb83e432_2
  - libzlib=1.3.2=h25fd6f3_3
  - ncurses=6.6=hdb14827_1
  - openjdk=17.0.18=h1602c4f_18
  - openssl=3.6.4=h781a0a9_0
  - pcre2=10.47=h8b3dc9c_1
  - pixman=0.46.4=h54a6638_3
  - pthread-stubs=0.4=hb03c661_1003
  - samtools=1.24=h9dcdb79_1
  - xorg-libice=1.1.2=h280c20c_0
  - xorg-libsm=1.2.6=h0d788c3_1
  - xorg-libx11=1.8.13=he1eb515_1
  - xorg-libxau=1.0.12=hb03c661_2
  - xorg-libxdmcp=1.1.5=hb03c661_2
  - xorg-libxext=1.3.7=h7cc23a3_1
  - xorg-libxfixes=6.0.2=h7cc23a3_1
  - xorg-libxi=1.8.3=h7cc23a3_1
  - xorg-libxrandr=1.5.5=h7cc23a3_1
  - xorg-libxrender=0.9.12=hb03c661_1
  - xorg-libxt=1.3.1=h7cc23a3_1
  - xorg-libxtst=1.2.5=h7cc23a3_4
  - zstd=1.5.7=hb78ec9c_7
```

### snakemake_env

Save the YAML block as environment.snakemake_env.yml, then run:

```bash
conda env create -f environment.snakemake_env.yml
```

```yaml
name: snakemake_env
channels:
  - conda-forge
  - bioconda
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - aioeasywebdav=2.4.0=pyha770c72_0
  - aiohappyeyeballs=2.7.1=pyhd8ed1ab_0
  - aiohttp=3.14.3=py310h31b6992_0
  - aiosignal=1.4.0=pyhd8ed1ab_0
  - alsa-lib=1.2.16.1=h7cc23a3_1
  - amply=0.1.7=pyhd8ed1ab_0
  - annotated-doc=0.0.5=pyhcf101f3_0
  - annotated-types=0.8.0=pyhd8ed1ab_0
  - appdirs=1.4.4=pyhd8ed1ab_1
  - argparse-dataclass=2.0.0=pyhd8ed1ab_1
  - async-timeout=5.0.1=pyhcf101f3_2
  - attrs=26.1.0=pyhcf101f3_0
  - backports.zstd=1.7.0=py310hff74037_0
  - bcrypt=5.0.0=py310hd8f68c5_1
  - boto3=1.43.83=pyhd8ed1ab_0
  - botocore=1.43.83=pyhd8ed1ab_0
  - brotli=1.2.0=h505cf86_3
  - brotli-bin=1.2.0=h9908984_3
  - brotli-python=1.2.0=py310h25d4438_3
  - bzip2=1.0.8=hda65f42_10
  - c-ares=1.34.8=hebe6cf0_2
  - ca-certificates=2026.7.22=hbd8a1cb_0
  - cairo=1.18.4=h3c89d7e_3
  - certifi=2026.7.22=pyhd8ed1ab_0
  - cffi=2.1.1=py310hd4ea0ea_2
  - charset-normalizer=3.5.1=pyhd8ed1ab_0
  - coin-or-cbc=2.10.13=h4d16d09_1
  - coin-or-cgl=0.60.10=hc46dffc_1
  - coin-or-clp=1.17.11=hc03379b_1
  - coin-or-osi=0.108.12=hf4fecb4_1
  - coin-or-utils=2.11.13=hc93afbd_1
  - coincbc=2.10.13=1_metapackage
  - colorama=0.4.6=pyhd8ed1ab_1
  - coloredlogs=15.0.1=pyhd8ed1ab_4
  - configargparse=1.7.5=pyhcf101f3_0
  - connection_pool=0.0.3=pyhd3deb0d_0
  - contourpy=1.3.2=py310h3788b33_0
  - cryptography=50.0.0=py310h20b05aa_1
  - cycler=0.12.1=pyhcf101f3_2
  - cyrus-sasl=2.1.28=hac629b4_1
  - datrie=0.8.3=py310h7c4b9e2_0
  - dbus=1.16.2=he8c428d_2
  - defusedxml=0.7.1=pyhd8ed1ab_0
  - docutils=0.23=pyhcf101f3_0
  - double-conversion=3.4.0=hecca717_0
  - dpath=2.2.0=pyha770c72_1
  - dropbox=12.0.2=pyhd8ed1ab_1
  - eido=0.2.6=pyhd8ed1ab_0
  - filechunkio=1.8=py_2
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.18.3=h4db4eae_1
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - fonttools=4.63.0=py310h3406613_0
  - formulaic=1.2.2=pyhd8ed1ab_0
  - freetype=2.14.3=ha770c72_2
  - frozenlist=1.8.0=py310h9548a50_0
  - ftputil=5.1.0=pyhd8ed1ab_0
  - gitdb=4.0.12=pyhd8ed1ab_0
  - gitpython=3.1.61=pyhcf101f3_0
  - google-api-core=2.34.0=pyhcf101f3_0
  - google-api-core-grpc=2.34.0=pyh17820f2_0
  - google-api-python-client=2.199.0=pyh332efcf_0
  - google-auth=2.57.0=pyhcf101f3_0
  - google-auth-httplib2=0.4.2=pyhcf101f3_0
  - google-cloud-core=2.7.0=pyhcf101f3_0
  - google-cloud-storage=3.13.1=pyhcf101f3_0
  - google-crc32c=1.8.0=py310hf432777_1
  - google-resumable-media=2.10.2=pyhcf101f3_0
  - googleapis-common-protos=1.75.2=pyhcf101f3_0
  - googleapis-common-protos-grpc=1.75.2=pyhcf101f3_0
  - graphite2=1.3.15=h54a6638_1
  - grpc-google-iam-v1=0.14.5=pyhcf101f3_0
  - grpcio=1.83.0=py310hc274405_0
  - grpcio-status=1.83.0=pyhcf101f3_0
  - h2=4.4.1=pyhcf101f3_0
  - hpack=4.2.0=pyhd8ed1ab_0
  - httplib2=0.32.0=pyhcf101f3_0
  - humanfriendly=10.0=pyh707e725_8
  - hyperframe=6.1.0=pyhd8ed1ab_0
  - icu=78.3=py310h44b86e0_2
  - idna=3.19=pyhcf101f3_0
  - interface_meta=1.3.0=pyhd8ed1ab_1
  - invoke=3.0.3=pyhd8ed1ab_0
  - jinja2=3.1.6=pyhcf101f3_1
  - jmespath=1.1.0=pyhcf101f3_1
  - jsonschema=4.26.0=pyhcf101f3_1
  - jsonschema-specifications=2025.9.1=pyhcf101f3_0
  - jupyter_core=5.9.1=pyhc90fa1f_0
  - keyutils=1.6.3=h7cc23a3_1
  - kiwisolver=1.5.1=py310he417ec5_0
  - krb5=1.22.2=hbc21106_2
  - lcms2=2.19.1=h9073bf1_2
  - ld_impl_linux-64=2.46.1=default_hbd61a6d_102
  - lerc=4.2.0=hdb68285_0
  - libabseil=20260526.0=cxx17_h0dc7533_2
  - libblas=3.11.0=10_h4a7cf45_openblas
  - libbrotlicommon=1.2.0=h39a168f_3
  - libbrotlidec=1.2.0=ha411449_3
  - libbrotlienc=1.2.0=h018ffa1_3
  - libcblas=3.11.0=10_h0358290_openblas
  - libclang-cpp23.1=23.1.0=default_h0acdd01_1
  - libclang13=23.1.0=default_h7855034_1
  - libcrc32c=1.1.2=h9c3ff4c_0
  - libcups=2.3.3=h7a8fb5f_6
  - libdeflate=1.25=hd45a770_1
  - libdrm=2.4.129=h7cc23a3_0
  - libedit=3.1.20250104=pl5321h373387f_1
  - libegl=1.7.0=ha4b6fd6_5
  - libegl-devel=1.7.0=ha4b6fd6_5
  - libexpat=2.8.1=hecca717_1
  - libffi=3.7.0=h81df57d_1
  - libfreetype=2.14.3=ha770c72_2
  - libfreetype6=2.14.3=h5e6c136_2
  - libgcc=16.2.0=ha9f2e26_4
  - libgcc-ng=16.2.0=h69a702a_4
  - libgfortran=16.2.0=h69a702a_4
  - libgfortran5=16.2.0=h6b99dfc_4
  - libgl=1.7.0=ha4b6fd6_5
  - libgl-devel=1.7.0=ha4b6fd6_5
  - libglib=2.88.3=he503a2a_2
  - libglvnd=1.7.0=ha4b6fd6_5
  - libglx=1.7.0=ha4b6fd6_5
  - libglx-devel=1.7.0=ha4b6fd6_5
  - libgomp=16.2.0=he0feb66_4
  - libgrpc=1.83.0=h792040b_0
  - libharfbuzz=14.4.0=h23af247_0
  - libiconv=1.18=h0cb94f2_3
  - libjpeg-turbo=3.2.0=hb03c661_1
  - liblapack=3.11.0=10_h47877c9_openblas
  - liblapacke=3.11.0=10_h6ae95b6_openblas
  - libllvm23=23.1.0=h474f4eb_0
  - liblzma=5.8.3=hb03c661_1
  - libnsl=2.0.1=hb9d3cd8_1
  - libntlm=1.8=hb9d3cd8_0
  - libopenblas=0.3.34=pthreads_hcf972fe_1
  - libopengl=1.7.0=ha4b6fd6_5
  - libpciaccess=0.19=hb03c661_1
  - libpng=1.6.58=h922cc85_1
  - libpq=18.6=h9d76c99_0
  - libprotobuf=7.35.1=h622638d_3
  - libre2-11=2025.11.05=h60473fc_2
  - libsodium=1.0.22=hebe6cf0_2
  - libsqlite=3.53.4=h13e7031_1
  - libstdcxx=16.2.0=h934c35e_4
  - libstdcxx-ng=16.2.0=hdf11a46_4
  - libtiff=4.7.2=hcc2c06a_1
  - libuuid=2.42.2=h5347b49_0
  - libvulkan-loader=1.4.357.0=h0e34353_2
  - libwebp-base=1.6.0=hd42ef1d_1
  - libxcb=1.17.0=hb83e432_2
  - libxcrypt=4.4.38=h280c20c_0
  - libxkbcommon=1.13.2=h51789e4_1
  - libxml2=2.15.3=h49c6c72_1
  - libxml2-16=2.15.3=hca6bf5a_1
  - libxslt=1.1.43=h711ed8c_1
  - libzlib=1.3.2=h25fd6f3_3
  - logmuse=0.3.0=pyhcf101f3_0
  - markdown-it-py=4.2.0=pyhd8ed1ab_0
  - markupsafe=3.0.3=py310h3406613_1
  - matplotlib=3.10.9=py310hff52083_0
  - matplotlib-base=3.10.9=py310hfde16b3_0
  - mdurl=0.1.2=pyhd8ed1ab_1
  - multidict=6.7.1=py310h3406613_0
  - munkres=1.1.4=pyhd8ed1ab_1
  - narwhals=2.25.0=pyhcf101f3_0
  - nbformat=5.11.1=pyhcf101f3_0
  - ncurses=6.6=hdb14827_1
  - numpy=2.2.6=py310hefbff90_0
  - oauth2client=4.1.3=pyhd8ed1ab_1
  - openjpeg=2.5.4=h55e1f5c_1
  - openldap=2.6.13=hbde042b_0
  - openssl=3.6.4=h781a0a9_0
  - packaging=26.3=pyhc364b38_0
  - pandas=2.3.3=py310h0158d43_2
  - paramiko=5.0.0=pyhd8ed1ab_0
  - patsy=1.0.3=pyhcf101f3_0
  - pcre2=10.47=h8b3dc9c_1
  - pephubclient=0.6.0=pyhd8ed1ab_0
  - peppy=0.40.8=pyhd8ed1ab_0
  - pillow=12.3.0=py310h5a73078_0
  - pip=26.2.1=pyh8b19718_0
  - pixman=0.46.4=h54a6638_3
  - platformdirs=4.11.5=pyhcf101f3_0
  - ply=3.11=pyhd8ed1ab_3
  - prettytable=3.18.0=pyhd8ed1ab_0
  - propcache=0.5.2=py310h3406613_0
  - proto-plus=1.28.4=pyhcf101f3_0
  - protobuf=7.35.1=py310hffa8988_1
  - psutil=7.2.2=py310h83e8816_2
  - pthread-stubs=0.4=hb03c661_1003
  - pulp=2.7.0=py310hff52083_1
  - pyasn1=0.6.4=pyhcf101f3_0
  - pyasn1-modules=0.4.2=pyhd8ed1ab_0
  - pycparser=3.0=pyhcf101f3_0
  - pydantic=2.13.5=pyhcf101f3_0
  - pydantic-core=2.46.5=py310hbfe5df0_1
  - pygments=2.21.0=pyhcf101f3_0
  - pynacl=1.6.2=py310h4e0494e_4
  - pyparsing=3.3.2=pyhcf101f3_0
  - pysftp=0.2.9=py_1
  - pyside6=6.11.2=py310ha8dd31c_0
  - pysocks=1.7.1=pyha55dd90_7
  - python=3.10.21=h267e890_0_cpython
  - python-dateutil=2.9.0.post0=pyhe01879c_2
  - python-fastjsonschema=2.22.2=pyhcf101f3_0
  - python-irodsclient=3.1.1=pyhff2d567_1
  - python-tzdata=2026.3=pyhd8ed1ab_0
  - python_abi=3.10=9_cp310
  - pytz=2026.3.post1=pyhcf101f3_0
  - pyu2f=0.1.5=pyhd8ed1ab_1
  - pyyaml=6.0.3=py310h3406613_1
  - qhull=2020.2=h434a139_5
  - qt6-main=6.11.2=pl5321h9df5c37_0
  - re2=2025.11.05=h94463f1_2
  - readline=8.3=hd6e31c0_1
  - referencing=0.37.0=pyhcf101f3_0
  - requests=2.34.2=pyhcf101f3_0
  - reretry=0.11.8=pyhd8ed1ab_1
  - rich=15.0.0=pyhcf101f3_0
  - rpds-py=0.30.0=py310hd8f68c5_0
  - rsa=4.9.1=pyhd8ed1ab_0
  - s3transfer=0.19.2=pyhd8ed1ab_0
  - scipy=1.15.2=py310h1d65ade_0
  - seaborn=0.13.2=hd8ed1ab_3
  - seaborn-base=0.13.2=pyhd8ed1ab_3
  - setuptools=81.0.0=pyh332efcf_0
  - setuptools-scm=10.2.1=pyhcf101f3_0
  - shellingham=1.5.4=pyhd8ed1ab_2
  - six=1.17.0=pyhe01879c_1
  - slacker=0.14.0=pyhd8ed1ab_1
  - smart_open=8.0.1=pyhcf101f3_0
  - smmap=5.0.3=pyhcf101f3_1
  - snakemake=7.32.4=hdfd78af_1
  - snakemake-minimal=7.32.4=pyhdfd78af_1
  - statsmodels=0.15.0=np2py310h83c7da3_1
  - stone=3.3.2=pyhd8ed1ab_0
  - stopit=1.1.2=pyhcf101f3_2
  - tabulate=0.10.0=pyhcf101f3_0
  - throttler=1.2.2=pyhd8ed1ab_0
  - tk=8.6.13=noxft_h1df4ec4_4
  - toml=0.10.2=pyhcf101f3_3
  - tomli=2.4.1=pyhcf101f3_0
  - toposort=1.10=pyhd8ed1ab_1
  - tornado=6.5.8=py310h7c4b9e2_0
  - traitlets=5.16.1=pyhcf101f3_0
  - typer=0.27.2=pyhcf101f3_0
  - typing-extensions=4.16.0=h69aa097_0
  - typing-inspection=0.4.4=pyhcf101f3_0
  - typing_extensions=4.16.0=pyhcf101f3_0
  - tzdata=2026c=h151e31d_0
  - ubiquerg=0.9.3=pyhd8ed1ab_0
  - unicodedata2=17.0.1=py310h7c4b9e2_0
  - uritemplate=4.2.0=pyhd8ed1ab_0
  - urllib3=2.7.0=pyhd8ed1ab_0
  - vcs_versioning=2.3.1=pyhcf101f3_0
  - wayland=1.26.0=hc1c935e_2
  - wcwidth=0.8.3=pyhcf101f3_0
  - wheel=0.48.0=pyhd8ed1ab_0
  - wrapt=2.4.0=py310h83e8816_0
  - xcb-util=0.4.1=h4f16b4b_2
  - xcb-util-cursor=0.1.6=hb03c661_0
  - xcb-util-image=0.4.0=hb711507_2
  - xcb-util-keysyms=0.4.1=hb711507_0
  - xcb-util-renderutil=0.3.10=hb711507_0
  - xcb-util-wm=0.4.2=hb711507_0
  - xkeyboard-config=2.48=h280c20c_0
  - xorg-libice=1.1.2=h280c20c_0
  - xorg-libsm=1.2.6=h0d788c3_1
  - xorg-libx11=1.8.13=he1eb515_1
  - xorg-libxau=1.0.12=hb03c661_2
  - xorg-libxcomposite=0.4.7=hb03c661_0
  - xorg-libxcursor=1.2.3=hb9d3cd8_0
  - xorg-libxdamage=1.1.6=hb9d3cd8_0
  - xorg-libxdmcp=1.1.5=hb03c661_2
  - xorg-libxext=1.3.7=h7cc23a3_1
  - xorg-libxfixes=6.0.2=h7cc23a3_1
  - xorg-libxi=1.8.3=h7cc23a3_1
  - xorg-libxrandr=1.5.5=h7cc23a3_1
  - xorg-libxrender=0.9.12=hb03c661_1
  - xorg-libxtst=1.2.5=h7cc23a3_4
  - xorg-libxxf86vm=1.1.7=hb03c661_0
  - xorg-xorgproto=2025.1=h280c20c_1
  - yaml=0.2.5=hebe6cf0_3
  - yarl=1.24.5=py310h3406613_0
  - yte=1.9.4=pyhd8ed1ab_0
  - zlib-ng=2.3.3=hce19668_1
  - zstd=1.5.7=hb78ec9c_7
```

