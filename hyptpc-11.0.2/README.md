k18geant4
=========

K1.8 geant4 simulation tool.

**Note that this README needs to be modified.**

## Platform

This tool is developed on the platform of KEKCC, RHEL 9.8.
- g++ (GCC) 11.5.0
- ROOT 6.32/04
- Geant4 11.2.2


## Anaconda setting

To use Python,
it is necessary to build the Anaconda local environment once using the `conda` command as follows.
Note that it is recommended to use `conda install` instead of `pip install` in the anaconda environment.

```sh
$ conda create -n myenv python=3.9 # myenv is an example name
$ conda activate myenv
$ conda install numpy psutil pyyaml rich
```

Add the following line in .bashrc to activate your environment.

```sh
conda activate myenv
```

If the prompt header of conda is annoying, add the following line in .condarc.

```yaml
changeps1: False
```


## How to install

Set environment variables.

```shell
. /sw/packages/root/6.32.04/bin/thisroot.sh
. /sw/packages/geant4/11.2.2/bin/geant4.sh
. /sw/packages/geant4/11.2.2/share/Geant4/geant4make/geant4make.sh
export MAKEFLAGS=-j40
conda activate myenv
```

then

CMakeLists.txt is removed from git management, so copy it from CMakeLists.txt.org first.

```shell
git clone --branch e72 git@github.com:hyptpc/g4-hyptpc.git
cd g4-hyptpc/hyptpc-11.0.2
cp CMakeLists.txt.org CMakeLists.txt
./build.sh
```



## How to use

Arguments of ConfFile and OutputName are necessary.
G4Macro is an optional argument.

```shell
./bin/G4HypTPC [ConfFile] [OutputName] (G4Macro)
./bin/G4HypTPC param/conf/e72_beam.conf foo.root
./bin/G4HypTPC param/conf/e72_beam.conf foo.root bar.mac
```



## Parameters

Some parameter files that are out of the git control should be linked.

```shell
ln -s /group/had/sks/E72/software/param/BEAM/* param/BEAM/
ln -s /group/had/sks/E72/software/fieldmap .
```
