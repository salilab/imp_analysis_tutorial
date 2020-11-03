Folder containing ihm deposition script and output mmcif file for IMP analysis tutorial 2020

### Files:
`create_ihm_cif_file.py` : Python script to build mmcif file.  Requires IMP. 
`tutorial_output.cif`    : Output of `create_ihm_cif_file.py`, an mmcif file that archives all four stages of modeling and output models

### Usage:

`python create_ihm_cif_file.py`

This will output the mmcif file `tutorial_output.cif`, along with `xl.db` files, which can be ignored. 

This .cif file is appropriate for deposition in the [PDB-Dev repository](https://pdb-dev.wwpdb.org/) for integrative models.

### Further info

A [complete tutorial](https://integrativemodeling.org/tutorials/deposition/) on archiving and depositing integrative models using IMP is available. The code can be [downloaded here](https://github.com/salilab/imp_deposition_tutorial).
