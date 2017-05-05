import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.proteomics
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros

import os


output_objects = [] # this will contain object that produce output log

# setting up topology

m = IMP.Model()

pdbfile="../../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb"
fasta_file="../../../Basteplate.crosslink.fasta"


conversion_chainid_name = {"A":"TssK"}

conversion_fastanames_name = {"TssKHis6":"TssK"}
                           
inv_conversion_chainid_name = {v: k for k, v in conversion_chainid_name.items()}
inv_conversion_fastanames_name = {v: k for k, v in conversion_fastanames_name.items()}

mol_list_2=["TssK"]

s = IMP.pmi.topology.System(m)
st = s.create_state()
seqs = IMP.pmi.topology.Sequences(fasta_file,conversion_fastanames_name)

pdbseqs=IMP.pmi.topology.PDBSequences(m,pdbfile,conversion_chainid_name)

#offsets=IMP.pmi.topology.fasta_pdb_alignments(seqs,pdbseqs,show=True)

offsets={"TssK":{(1,451):0}}

mols={}
for n,name in enumerate(mol_list_2):
    mol = st.create_molecule(                      # create molecule
        name,
        sequence=seqs[name],
        chain_id=inv_conversion_chainid_name[name])
    atomic_set=IMP.pmi.tools.OrderedSet()
    for group in offsets[name]:
        print group, offsets[name][group]
        atomic = mol.add_structure(pdbfile,
                                chain_id=inv_conversion_chainid_name[name],
                                res_range=group,
                                offset=offsets[name][group],soft_check=True)
        mol.add_representation(atomic,                 # res 1,10 for structured regions
                               resolutions=[1],
                               color=0.0)
        atomic_set|=atomic
    if len(mol[:]-atomic_set) > 0:
        mol.add_representation(mol[:]-atomic_set,          # res 10 for unstructured regions
                               resolutions=[10],
                               color=0.0)
    mols[name]=mol

# calling System.build() creates all States and Molecules (and their representations)
#  Once you call build(), anything without representation is destroyed.
#  You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
#  However these functions will only return BUILT representations
root_hier = s.build()

# Setup degrees of freedom
#  The DOF functions automatically select all resolutions
#  Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(m)



dof.create_rigid_body(mols.values(),
                      max_trans=0.0,
                      max_rot=0.0,
                      nonrigid_max_trans=4.0)

# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)

mols = [mols[mol] for mol in mols]
crs = []
for mol in mols:
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    output_objects.append(cr)
    crs.append(cr)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.add_to_model()
output_objects.append(evr)

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    #crosslink_restraints=xls,                     # will display like XLs
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='output_TssK/',
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=100)                   # increase number of frames to get better results!
rex.execute_macro()

