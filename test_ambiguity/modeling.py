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
import IMP.pmi.restraints.em
import IMP.pmi.plotting
import IMP.pmi.plotting.topology
import IMP.bayesianem
import IMP.bayesianem.restraint
import tempfile,os
import sys
import IMP.em
import numpy as np
import math 


output_objects=[]
mdl = IMP.Model()

'''
#Sequence alignment
seqs = IMP.pmi.topology.Sequences("../data/Basteplate.crosslink.fasta")

pdbseqs=IMP.pmi.topology.PDBSequences(mdl,"../data/pdb/TssK_premilinar_xray.pdb")

offsets=IMP.pmi.topology.fasta_pdb_alignments(seqs,pdbseqs,show=True)
'''

###################### SYSTEM SETUP #####################
# Read sequences etc

topology='''
|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|
|TssF   |blue       |../data/Basteplate.crosslink.fasta|StrepTssF|BEADS||1,40||20|10||||
|TssF.1 |cyan       |../data/Basteplate.crosslink.fasta|StrepTssF|BEADS||1,40||20|10||||
'''


# Normally the topology table is kept in a text file but here we just write it to a temporary one
tf = tempfile.NamedTemporaryFile(delete=False)
tf.write(topology)
tf.close()

# The TopologyReader reads the text file, and the BuildSystem macro constructs it

reader = IMP.pmi.topology.TopologyReader(tf.name,
                                         pdb_dir = '.',
                                         fasta_dir = '.',
                                         gmm_dir = 'gmm')
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(reader) # note you can call this multiple times to create a multi-state system

#IMP.pmi.plotting.topology.draw_component_composition(bs)


hier, dof = bs.execute_macro(max_rb_trans=2.0, max_rb_rot=0.1, max_bead_trans=2.0, max_srb_trans=2.0,max_srb_rot=0.1)

#particles=IMP.atom.Selection(hier,molecule="TssF",residue_indexes=range(1,315+1)).get_selected_particles()
#dof.disable_movers(particles)

#IMP.core.XYZ(particles[0]).set_coordinates((0,0,0))
#IMP.core.XYZ(particles[1]).set_coordinates((40,0,0))

##parts_sub2=IMP.atom.Selection(hier,molecule="Sub2").get_selected_particles()
##parts_yra1=IMP.atom.Selection(hier,molecule="Yra1",residue_indexes=range(208,225)).get_selected_particles()

##fixed_xyz,fixed_rb=dof.fix_particles(parts_sub2+parts_yra1,IMP.core.RigidBodyMover)


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
for molname in moldict:
    for mol in moldict[molname]:
        IMP.pmi.tools.display_bonds(mol)
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=2.0)
        cr.add_to_model()
        cr.set_label(molname)
        output_objects.append(cr)
        crs.append(cr)
        mols.append(mol)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.add_to_model()
evr.set_weight(0.5)
output_objects.append(evr)



import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
from IMP.pmi.io.crosslink import FilterOperator as FO
import IMP.crosslinkms
import IMP.crosslinkms.restraint_pmi2
import operator


rplp=IMP.pmi.io.crosslink.ResiduePairListParser("MSSTUDIO")
cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter(rplp)
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_site_pairs_key("Selected Sites")
cldbkc.set_id_score_key("E")
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.set_name("master")

for csv in ["test.csv"]:

#for csv in ["../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFG.csv","../data/crosslinks/old/Exp2-KFG-081016-2.csv"]:

    cldb_tmp=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    
    cldb_tmp.create_set_from_file(csv)
    
    cldb_tmp.set_name(csv)
    cldb.append_database(cldb_tmp)
    
    
for xl in cldb:
    xl[cldb.id_score_key]=-math.log(xl[cldb.id_score_key])

cldb.create_new_keyword(cldb.sigma1_key)
cldb.create_new_keyword(cldb.sigma2_key)

cldb.set_value(cldb.sigma1_key,'SIGMA')
cldb.set_value(cldb.sigma2_key,'SIGMA')
#cldb.set_value(cldb.sigma1_key,'FLEXIBLE-1',FO(cldb.residue1_links_number_key,operator.gt,5))
#cldb.set_value(cldb.sigma2_key,'FLEXIBLE-1',FO(cldb.residue2_links_number_key,operator.gt,5))
#cldb.set_value(cldb.sigma1_key,'FLEXIBLE-2',FO(cldb.residue1_links_number_key,operator.gt,30))
#cldb.set_value(cldb.sigma2_key,'FLEXIBLE-2',FO(cldb.residue2_links_number_key,operator.gt,30))


#for xl in cldb:
#    print xl[cldb.sigma1_key],xl[cldb.residue1_links_number_key],xl[cldb.sigma2_key],xl[cldb.residue2_links_number_key]


data=IMP.crosslinkms.DataFixedLength("/Users/riccardopellarin/imp-project/imp-060516/imp/modules/crosslinkms/data/data.21.10000.gaussian.1.rmf")
xlmovers=[]
xls=[]

xl = IMP.crosslinkms.restraint_pmi2.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                CrossLinkDataBase=cldb,
                                length=21.0,
                                slope=0.02,
                                resolution=1.0,
                                label="XL",
                                data=data,
                                alpha=10,
                                beta=100,
                                isotonic_regression=False,#,cldb.id_score_key,
                                outlier=True) #cldb.id_score_key)
xl.add_to_model()

sigma=xl.sigma_dictionary["SIGMA"][0]
sigma.set_scale(11.0)
xl.set_sigma_is_sampled(True)
#sigma=xl.sigma_dictionary["FLEXIBLE-1"][0]
#sigma.set_scale(11.0)
#sigma=xl.sigma_dictionary["FLEXIBLE-2"][0]
#sigma.set_scale(11.0)
#xl.set_sigma_is_sampled(True)
psi=xl.psi_dictionary["PSI"][0]
psi.set_scale(0.05)
xl.set_psi_is_sampled(True)
output_objects.append(xl)
xls.append(xl)
xlmovers+=xl.get_movers()




'''
IMP.pmi.tools.shuffle_configuration(hier,
                                        max_translation=100,avoidcollision_rb=True, avoidcollision_fb=True)


'''

output_dir='./'
output_prefix=''
## Run replica exchange Monte Carlo sampling
rex_obj=None
weights = 1.0 * np.array([1, 1]*100)
i=0;
for w in weights:
    #gemt.set_weight(w)
    print rex_obj
    num_frames=100
    if w==1:
        num_frames=100000
    if rex_obj is None:
        rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                            root_hier=hier,                          # pass the root hierarchy
                                            crosslink_restraints=xls,                     # will display like XLs
                                            monte_carlo_sample_objects=dof.get_movers()+xlmovers,# pass MC movers
                                            global_output_directory='%s/%soutput_%d_%g/'%(output_dir,output_prefix, i,w),
                                            output_objects=output_objects,
                                            monte_carlo_steps=10,
                                            self_adaptive=True,
                                            number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                            number_of_frames=num_frames)                   # increase number of frames to get better results!
        rex.execute_macro()
        rex_obj=rex.replica_exchange_object
    else:
        rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                                root_hier=hier,                          # pass the root hierarchy
                                                crosslink_restraints=xls,                     # will display like XLs
                                                monte_carlo_sample_objects=dof.get_movers()+xlmovers,# pass MC movers
                                                global_output_directory='%s/%soutput_%d_%g/'%(output_dir,output_prefix,i,w),
                                                output_objects=output_objects,
                                                monte_carlo_steps=10,
                                                number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                                number_of_frames=num_frames,                   # increase number of frames to get better results!
                                                replica_exchange_object=rex_obj)
        rex.execute_macro()
    i=i+1