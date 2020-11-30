### Installation of CPSR using Conda

This is an alternative installation approach that does not require Docker on your machine. At the moment it works only for Linux systems.

A prerequisite is that you have Conda installed. First download the Conda package manager. Get it with:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
```

Run the following to add Conda into your PATH. You can even put that into your `~/.bashrc` or `~/.zshrc` to avoid re-running this in the future:

```
. ./miniconda/etc/profile.d/conda.sh
```

#### Alternative 1
Create a new environment (`-n cpsr`), install the _cpsr_ Conda package directly from Anaconda Cloud:

```
conda create -n cpsr -c conda-forge -c bioconda -c pcgr cpsr

```

#### Alternative 2
Build the _cpsr_ package from source, which is useful if you need to use the development code from the repository:

```
conda install conda-build
export CHANNELS="-c conda-forge -c bioconda -c pcgr -c defaults"
conda build $CHANNELS conda_pkg/cpsr
conda install --use-local $CHANNELS cpsr
```

For both alternatives you also need to download the reference data bundle for your genome build:

```
wget http://insilico.hpc.uio.no/pcgr/pcgr.databundle.20201123.grch37.tar.gz -O grch37.tar.gz
wget http://insilico.hpc.uio.no/pcgr/pcgr.databundle.20201123.grch38.tar.gz -O grch38.tar.gz
tar -xzf grch37.tar.gz  # will extract into ./data/grch37/
tar -xzf grch38.tar.gz  # will extract into ./data/grch38/
```

There is a chance you'll encounter errors during the installation. Due to ongoing updates of the packages in public repositories, some packages might end up conflicting with each other or missing for your system. So try to stick to the dockerized version of CPSR whenever possible.

### Running condarized CPSR

Activate your environment with:

```
conda activate cpsr
```

Run CPSR with `--no-docker` flag. The `--pcgr_dir` argument now doesn't have to contain anything but a `data` directory that you downloaded.

If you encounter errors with VEP, you may need to unset/reset PERL5LIB (e.g. `export PERL5LIB=""`), see the [following issue](https://github.com/bioconda/bioconda-recipes/issues/4390)
