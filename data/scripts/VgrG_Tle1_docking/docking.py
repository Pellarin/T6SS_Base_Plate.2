import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import tempfile,os
import sys


output_objects=[]
sample_objects=[]


###################### SYSTEM SETUP #####################
# Read sequences etc
topology='''
|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|
|VgrGA  |blue      |Basteplate.crosslink.fasta|VgrG-II|../../modeller/VgrG/VGRG.B99990001.pdb|A|1,841|+8|5|0 |2 |  | |
|VgrGB  |blue      |Basteplate.crosslink.fasta|VgrG-II|../../modeller/VgrG/VGRG.B99990001.pdb|B|842,1682|-833|5|0 |2 |  | |
|VgrGC |blue      |Basteplate.crosslink.fasta|VgrG-II|../../modeller/VgrG/VGRG.B99990001.pdb|C|1683,2523|-1674|5|0 |2 |  | |
|Tle1  |green     |Basteplate.crosslink.fasta|Tle1-His|c4o5pB_.1.pdb|A|1,END|0|5|0 |1 |  | |
'''


# Normally the topology table is kept in a text file but here we just write it to a temporary one
tf = tempfile.NamedTemporaryFile(delete=False)
tf.write(topology)
tf.close()

# The TopologyReader reads the text file, and the BuildSystem macro constructs it
mdl = IMP.Model()
reader = IMP.pmi.topology.TopologyReader(tf.name,
                                         pdb_dir = './',
                                         fasta_dir = '../../',
                                         gmm_dir = './')
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(reader) # note you can call this multiple times to create a multi-state system
hier, dof = bs.execute_macro(max_rb_rot=0.2)

# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
for molname in moldict:
    for mol in moldict[molname]:
        IMP.pmi.tools.display_bonds(mol)
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        cr.set_label(molname)
        output_objects.append(cr)
        crs.append(cr)
        mols.append(mol)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.add_to_model()
output_objects.append(evr)

import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
from IMP.pmi.io.crosslink import FilterOperator as FO
import operator


rplp=IMP.pmi.io.crosslink.ResiduePairListParser("MSSTUDIO")
cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter(rplp)
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_site_pairs_key("Selected Sites")
cldbkc.set_id_score_key("E")
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.set_name("master")

for csv in ["../../crosslinks/Vgrg-Tle1.csv","../../crosslinks/Vgrg-Tle1-2-1008.csv"]:
    cldb_tmp=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb_tmp.create_set_from_file(csv)
    cldb_tmp.clone_protein("VgrG-II","VgrGB")
    cldb_tmp.clone_protein("VgrG-II","VgrGC")
    cldb_tmp.rename_proteins({"VgrG-II":"VgrGA"})
    cldb_tmp.rename_proteins({"Tle1-His":"Tle1"})
    cldb_tmp.set_name(csv)
    cldb.append_database(cldb_tmp)

cldb.create_new_keyword("Redundancy_for_filter",values_from_keyword="Redundancy")
fo=FO("Decision",operator.eq,"Accepted")&FO("Redundancy_for_filter",operator.gt,2)
fcldb=cldb.filter(fo)


xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                            CrossLinkDataBase=fcldb,
                            length=21.0,
                            slope=0.01,
                            resolution=1.0,
                            label="XL")
xl.add_to_model()
sigma=xl.sigma_dictionary["SIGMA"][0]
psi=xl.psi_dictionary["PSI"][0]
sigma.set_scale(11.0)
psi.set_scale(0.05)
xl.set_psi_is_sampled(True)
xl.set_sigma_is_sampled(True)
sample_objects.append(xl)
output_objects.append(xl)



###################### SAMPLING #####################
# First shuffle the system
IMP.pmi.tools.shuffle_configuration(hier,
                                    max_translation=100)

# Quickly move all flexible beads into place
dof.optimize_flexible_beads(10)

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          # pass the root hierarchy
                                    crosslink_restraints=[xl],                     # will display like XLs
                                    monte_carlo_sample_objects=dof.get_movers()+xl.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=output_objects,
                                    monte_carlo_temperature=1.0,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    monte_carlo_steps=10,
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=10000)                   # increase number of frames to get better results!
rex.execute_macro()
