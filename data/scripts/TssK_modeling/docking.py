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
import IMP.bayesianem
import IMP.bayesianem.restraint
import tempfile,os
import sys
import math


output_objects=[]
sample_objects=[]


###################### SYSTEM SETUP #####################
# Read sequences etc
topology='''
|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|
|TssK  |blue   |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|1,184||5|10|1|1||
|TssK  |blue   |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|185,193||5|10||1||
|TssK  |blue   |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|194,315||5|10|2|1||
|TssK  |blue   |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|316,317||5|10||1||
|TssK  |blue   |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|318,451||5|10|3|1||
|TssK.1  |green        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|1,184||5|10||||
|TssK.1  |green        |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|185,193||5|10||||
|TssK.1  |green        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|194,315||5|10||||
|TssK.1  |green        |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|316,317||5|10||||
|TssK.1  |green        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|318,451||5|10||||
|TssK.2  |yellow        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|1,184||5|10||||
|TssK.2  |yellow        |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|185,193||5|10||||
|TssK.2  |yellow        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|194,315||5|10||||
|TssK.2  |yellow        |Basteplate.crosslink.fasta|TssKHis6|BEADS|A|316,317||5|10||||
|TssK.2  |yellow        |Basteplate.crosslink.fasta|TssKHis6|../../xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|318,451||5|10||||
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
        mols.append(mol)
        if IMP.atom.Copy(mol.hier).get_copy_index() == 0: 
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
            cr.add_to_model()
            cr.set_label(molname)
            output_objects.append(cr)
            crs.append(cr)
        

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.add_to_model()
evr.set_weight(0.1)
output_objects.append(evr)


copies=bs.get_molecules()[0]['TssK']

rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0)


for n,c in enumerate(copies[1:]):
    print(n,c)
    rotation_angle = 2.0 * math.pi / float(len(copies)) * float(n+1)
    rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
    transformation3D =IMP.algebra.get_rotation_about_point(IMP.algebra.Vector3D(126.41989953, 126.03652216,   0.00000000), rotation3D)

    dof.constrain_symmetry(
                       [copies[0]],
                       [c],
                       transformation3D,
                       resolution='all')


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

for csv in ["../../crosslinks/Exp1-BP-CXMS-12aout16-KFG.csv","../../crosslinks/Exp2-KFG-081016-2.csv",
            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGE.csv","../../crosslinks/Exp2-KFGE-081016-2.csv",
            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGV.csv","../../crosslinks/Exp3-KFG-081016-3.csv",
            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGVE.csv","../../crosslinks/Exp3-KFGE-081016-3.csv",
            "../../crosslinks/Vgrg-Tle1.csv","../../crosslinks/Vgrg-Tle1-2-1008.csv"]:

    cldb_tmp=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    
    cldb_tmp.create_set_from_file(csv)
    
    cldb_tmp.rename_proteins({"TssKHis6":"TssK"})
    cldb_tmp.set_name(csv)
    cldb.append_database(cldb_tmp)

fo=FO("Decision",operator.eq,"Accepted")

fcldb=cldb.filter(fo)


xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                            CrossLinkDataBase=fcldb,
                            length=21.0,
                            slope=0.02,
                            resolution=1.0,
                            label="XL")
xl.add_to_model()
sigma=xl.sigma_dictionary["SIGMA"][0]
psi=xl.psi_dictionary["PSI"][0]
sigma.set_scale(10.0)
psi.set_scale(0.05)
xl.set_psi_is_sampled(False)
xl.set_sigma_is_sampled(True)
sample_objects.append(xl)
output_objects.append(xl)

densities = IMP.atom.Selection(hier,molecule="TssK",representation_type=IMP.atom.DENSITIES).get_selected_particles()


gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                 '../TssK.50.txt',
                                                 target_mass=147477,
                                                 slope=0.01,
                                                 weight=0.5,
                                                 target_is_rigid_body=False)
gem.add_to_model()
gem.set_label("TssK")
output_objects.append(gem)


###################### SAMPLING #####################
# First shuffle the system
#IMP.pmi.tools.shuffle_configuration(hier,
#                                    max_translation=100)

# Quickly move all flexible beads into place
#dof.optimize_flexible_beads(10)

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          # pass the root hierarchy
                                    #crosslink_restraints=[xl],                     # will display like XLs
                                    monte_carlo_sample_objects=dof.get_movers(),#+xl.get_movers(),  # pass MC movers
                                    global_output_directory='output_test/',
                                    output_objects=output_objects,
                                    monte_carlo_temperature=1.0,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    simulated_annealing=True,
                                    monte_carlo_steps=10,
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=10000)                   # increase number of frames to get better results!
rex.execute_macro()
