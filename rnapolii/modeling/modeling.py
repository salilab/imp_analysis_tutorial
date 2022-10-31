"""
#############################################
##  IMP Tutorial Script
##
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: Papers where this data is shown...
#
"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import ihm.cross_linkers

import os
import sys


#---------------------------
# 1. Define Input Data and Output Directories
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"
output_directory = sys.argv[1]
output_index = sys.argv[2]

#--------------------------
# 2. Scoring Parameters
#--------------------------

#--------------
# ----- Sterochemistry and Physical Restraints
ev_weight = 1.0           # Weight of excluded volume restraint
connectivity_scale = 2.0  # weight of Connectivity restraint

#--------------
# ----- Electron Microscopy Restraint
em_weight = 80.0          # EM Restraint weight
em_slope = 0.000001       # EM restraint slope (to direct model into EM density)

# Gaussian Mixture Model (GMM) file created from EM map:
# python create_gmm.py emd_1883.map.mrc 50 emd_1883.gmm50.txt 
target_gmm_file = datadirectory+'emd_1883.gmm50.txt'

#--------------
# ----- Crosslinking Restraints
xl_slope = 0.002

# -----
# -- Crosslink set 1
xl1_file = datadirectory+'polii_xlinks.csv'

# Define the column headers in the XL1 file:
xldbkwc1 = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc1.set_protein1_key("pep1.accession")
xldbkwc1.set_protein2_key("pep2.accession")
xldbkwc1.set_residue1_key("pep1.xlinked_aa")
xldbkwc1.set_residue2_key("pep2.xlinked_aa")

xl1_linker = ihm.cross_linkers.dss # We're using a DSS crosslinker
xl1_length = 21.0                  # Restraint length is defined as 21.0 Angstroms
xl1_weight = 1.0                   # The weight of this restraint is 1.0 

# -----
# -- Crosslink Set 2
xl2_file = datadirectory+'polii_juri.csv'

# Define the column headers in the XL file:
xldbkwc2 = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc2.set_protein1_key("prot1")
xldbkwc2.set_protein2_key("prot2")
xldbkwc2.set_residue1_key("res1")
xldbkwc2.set_residue2_key("res2")

xl2_linker = ihm.cross_linkers.bs3 # BS3 crosslinker
xl2_length = 21.0
xl2_weight = 1.0

#--------------------
# 3. Sampling Parameters
#--------------------
num_frames = int(sys.argv[3])  # Number of frames in MC run
num_mc_steps = 10         # Number of MC steps per frame
mc_temperature = 1.0      # Temperature for MC

# --- Simulated Annealing (sa)
#  - Alternates between two MC temperatures
sim_annealing = True    # If true, run simulated annealing
sa_min_temp_steps = 100 # Steps at min temp
sa_max_temp_steps = 20  # Steps at max temp
sa_temps = (1.0, 2.5)   # Sim annealing temperatures

# Replica Exchange (rex)
rex_temps = (1.0, 2.5)  # Temperature bounds for replica exchange

#----------------
#------------------------------
# Here is where the real work begins...
#------------------------------
#----------------

# Initialize model
m = IMP.Model()

# Read in the topology file.
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=datadirectory,
                                  fasta_dir=datadirectory,
                                  gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Each state can be specified by a topology file.
bs.add_state(topology)

# Build the system representation and degrees of freedom
root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)

# Fix all rigid bodies but not Rpb4 and Rpb7 (the stalk)
# First select and gather all particles to fix.
fixed_particles=[]
for prot in ["Rpb1","Rpb2","Rpb3","Rpb5","Rpb6","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"]:
    fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot).get_selected_particles()

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,
                                          IMP.pmi.TransformMover])

# Randomize the initial configuration before sampling, of only the molecules
# we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = [] # reporter objects (for stat files)

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files

# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=connectivity_scale, label=molname)
    cr.add_to_model()
    outputobjects.append(cr)

# Excluded Volume Restraint
#  To speed up this expensive restraint, we evaluate it at resolution 10
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.set_weight(ev_weight)
ev.add_to_model()
outputobjects.append(ev)


# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("pep1.accession")
xldbkwc.set_protein2_key("pep2.accession")
xldbkwc.set_residue1_key("pep1.xlinked_aa")
xldbkwc.set_residue2_key("pep2.xlinked_aa")

xl1db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc1)
xl1db.create_set_from_file(xl1_file)

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl1db,
                                   length=xl1_length,
                                   slope=xl_slope,
                                   resolution=1.0, # Evaluate crosslink at resolution 1 (CA atom position)
                                   label="Trnka",
                                   linker=xl1_linker,
                                   weight=xl1_weight)

xl1.add_to_model()             # crosslink must be added to the model
outputobjects.append(xl1)


# Crosslinks - dataset 2
#  We can easily add a second set of crosslinks.
#  These have a different format and label, but other settings are the same

xl2db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc2)
xl2db.create_set_from_file(xl2_file)

xl2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl2db,
                                   length=xl2_length,
                                   slope=xl_slope,
                                   resolution=1.0,
                                   label="Chen",
                                   linker=xl2_linker,
                                   weight=xl2_weight)
xl2.add_to_model()
outputobjects.append(xl2)


# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-Bayesian
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                 target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=em_slope,
                                                 weight=em_weight)
gemt.add_to_model()
outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------

if '--test' in sys.argv: num_frames=100

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hier,
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    output_objects=outputobjects,
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=min(sa_temps),
                                    simulated_annealing_maximum_temperature=max(sa_temps),
                                    simulated_annealing_minimum_temperature_nframes=sa_min_temp_steps,
                                    simulated_annealing_maximum_temperature_nframes=sa_max_temp_steps,
                                    replica_exchange_minimum_temperature=min(rex_temps),
                                    replica_exchange_maximum_temperature=max(rex_temps),
                                    number_of_best_scoring_models=0,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory=output_directory+str(output_index))

# Start Sampling
mc1.execute_macro()
