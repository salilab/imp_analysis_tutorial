# Modeling of RNA Polymerase II stalk

To build the representation, scoring function and initiate the sampling protocol run:

`./run_rnapolii_modeling.sh out_dir N n_steps`

where `out_dir` is the prefix of the output directories. `N` is the number
of independent runs you wish to run (>2), and `n_steps` is the number of Monte Carlo steps to
run during one simulation. 

It is best to ensure that `N` is less than the number of cores on your machine.

Each 1000 steps will take on the order of a few minutes to run, depending on your machine.

running:

`time python modeling.py test 0 1000`

will give you a decent estimation of that time. 


