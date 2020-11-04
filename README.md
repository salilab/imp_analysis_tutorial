# imp_analysis_tutorial

This tutorial demonstrates the analysis of integrative modeling results using the [Integrative Modeling Platform](https://integrativemodeling.org) (IMP).

[./rnapolii](./rnapoliii) - Improved analysis routine for the [PMI modeling tutorial](https://integrativemodeling.org/tutorials/rnapolii_stalk/) of the RNA Polymerase III stalk.

Citations: In review

## Installation and Dependencies

### Getting the tutorial code and data
This tutorial can be downloaded from this page or the command line using `git clone github.com/salilab/imp_analysis_tutorial`

### Installing IMP and dependencies

The tutorial requires IMP v2.13 (April 2020) or later to run correctly.  

---Using anaconda---
The simplest way to install and run IMP and this tutorial is using [anaconda](https://www.anaconda.com/products/individual). 

Once installed, IMP and required dependencies can be installed using the following commands
```
conda install -c salilab imp pyrmsd
conda install numpy scikit-learn matplotlib pandas
conda install -c conda-forge hdbscan
```

---Using stand-alone IMP---
IMP [binaries or source code](https://integrativemodeling.org/download.html) are available for download. 

Python 3.X is recommended to run IMP.  Dependencies can be installed via `pip`:

```
pip install numpy scikit-learn matplotlib pandas hdbscan
```

A modified version of pyRMSD must be installed as well:
* [pyRMSD](https://github.com/salilab/pyRMSD)

## Running the tutorial

Navigate to `./rnapolii` and continue from the [rnapolii README](./rnapolii/README.md)
