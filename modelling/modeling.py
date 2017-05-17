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
|TssK   |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_A.pdb|A|1,315||20|10|1|||
|TssK   |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|2|||
|TssK.1 |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_A.pdb|A|1,315||20|10|3|||
|TssK.1 |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|4|||
|TssK.2 |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_A.pdb|A|1,315||20|10|5|||
|TssK.2 |red      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|6|||
|TssK.3 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_B.pdb|A|1,315||20|10|7|||
|TssK.3 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|8|||
|TssK.4 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_B.pdb|A|1,315||20|10|9|||
|TssK.4 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|10|||
|TssK.5 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/pdb/TssK_premilinar_xray_B.pdb|A|1,315||20|10|11|||
|TssK.5 |orange      |../data/Basteplate.crosslink.fasta|TssKHis6|../data/xraptor/TssK.197209.all_in_one/models/197209_model_1_chain.pdb|A|316,END||20|10|12|||
|TssF   |blue       |../data/Basteplate.crosslink.fasta|StrepTssF|BEADS||1,END||20|10|100|4,7,8|1|
|TssF.1   |cyan     |../data/Basteplate.crosslink.fasta|StrepTssF|BEADS||1,END||20|10|101|5,7,8| |
|TssG   |green       |../data/Basteplate.crosslink.fasta|TssGFlag|BEADS||1,END||20|10|102|6,8|3|
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


hier, dof = bs.execute_macro(max_rb_trans=4.0, max_rb_rot=0.3, max_bead_trans=4.0, max_srb_trans=4.0,max_srb_rot=0.3)



particles=IMP.atom.Selection(hier,molecule="TssK",residue_indexes=range(1,315+1)).get_selected_particles()
fixedxyz,fixedrbs=dof.disable_movers(particles,[IMP.core.RigidBodyMover])

for xyz in fixedxyz:
    IMP.core.XYZ(xyz).set_coordinates_are_optimized(False)
    
for rb in fixedrbs:
    IMP.core.RigidBody(rb).set_coordinates_are_optimized(False)


particles=IMP.atom.Selection(hier,molecule="TssF",copy_index=0).get_selected_particles()
for p in particles:
    rad=IMP.core.XYZR(p).get_radius()
    IMP.core.XYZR(p).set_radius(rad*0.75)
    p.set_value(IMP.FloatKey(4),296)
    p.set_value(IMP.FloatKey(5),216)
    p.set_value(IMP.FloatKey(6),262)

particles=IMP.atom.Selection(hier,molecule="TssF",copy_index=1).get_selected_particles()
for p in particles:
    rad=IMP.core.XYZR(p).get_radius()
    IMP.core.XYZR(p).set_radius(rad*0.75)
    p.set_value(IMP.FloatKey(4),278)
    p.set_value(IMP.FloatKey(5),233)
    p.set_value(IMP.FloatKey(6),182)

particles=IMP.atom.Selection(hier,molecule="TssG").get_selected_particles()
for p in particles:
    rad=IMP.core.XYZR(p).get_radius()
    IMP.core.XYZR(p).set_radius(rad*0.75)
    p.set_value(IMP.FloatKey(4),336)
    p.set_value(IMP.FloatKey(5),257)
    p.set_value(IMP.FloatKey(6),222)


##parts_sub2=IMP.atom.Selection(hier,molecule="Sub2").get_selected_particles()
##parts_yra1=IMP.atom.Selection(hier,molecule="Yra1",residue_indexes=range(208,225)).get_selected_particles()

##fixed_xyz,fixed_rb=dof.fix_particles(parts_sub2+parts_yra1,IMP.core.RigidBodyMover)


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
connected_pairs=[]
for molname in moldict:
    for n,mol in enumerate(moldict[molname]):        
        IMP.pmi.tools.display_bonds(mol)
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=1.0)
        cr.add_to_model()
        cr.set_label(molname)
        output_objects.append(cr)
        crs.append(cr)
        connected_pairs+=cr.get_particle_pairs()
        mols.append(mol)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = hier, resolution=10.0)
evr.add_to_model()
evr.set_weight(1.0)
evr.add_excluded_particle_pairs(connected_pairs)
output_objects.append(evr)




copies=bs.get_molecules()[0]['TssK']


'''
From Chimera
  Axis   0.98477941  -0.15168411   0.08486137
  Axis point   0.00000000 316.37920385 239.27469193
'''


rotational_axis=IMP.algebra.Vector3D(0.98477941,-0.15168411,0.08486137)
rotational_point=IMP.algebra.Vector3D(0.00000000,316.37920385,239.27469193)


copies=bs.get_molecules()[0]['TssK']


for n,c in enumerate([copies[1],copies[2]]):
    print(n,c)
    rotation_angle = 2.0 * math.pi / 3 * float(n+1)
    rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
    transformation3D =IMP.algebra.get_rotation_about_point(rotational_point, rotation3D)

    dof.constrain_symmetry(
                       [copies[0]],
                       [c],
                       transformation3D,
                       resolution='all')

'''
From Chimera
  Axis   -0.67968399  -0.53767865   0.49893020
  Axis point   0.00000000  40.27421766 426.23284010
'''

rotational_axis=IMP.algebra.Vector3D(-0.67968399,-0.53767865,0.49893020)
rotational_point=IMP.algebra.Vector3D(0.00000000,40.27421766,426.23284010)

for n,c in enumerate([copies[4],copies[5]]):
    print(n,c)
    rotation_angle = 2.0 * math.pi / 3 * float(n+1)
    rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
    transformation3D =IMP.algebra.get_rotation_about_point(rotational_point, rotation3D)

    dof.constrain_symmetry(
                       [copies[3]],
                       [c],
                       transformation3D,
                       resolution='all')






copies=bs.get_molecules()[0]['TssF']



#check this
#dof.constrain_symmetry(
#                       [copies[0]],
#                       [copies[1]],
#                       IMP.algebra.Transformation3D(IMP.algebra.Vector3D(100,0,0)), #IMP.algebra.get_identity_transformation_3d(),
#                       resolution='all',
#                       type="RIGID_BODY")


#setup the same internal coordinates




ps_constr=[]
lcp = IMP.container.ListSingletonContainer(mdl)
#smp = _InternalCoordinatesConstraint(mdl)
smp= IMP.pmi.InternalCoordinatesConstraint()
for (pref,pconstr) in zip(IMP.atom.get_leaves(copies[0].hier),IMP.atom.get_leaves(copies[1].hier)):
    pi_ref=pref.get_particle_index()
    pi_constr=pconstr.get_particle_index()
    ps_constr.append(pconstr.get_particle())
    
    if IMP.core.NonRigidMember.get_is_setup(mdl,pi_ref) and IMP.core.NonRigidMember.get_is_setup(mdl,pi_constr):
        
        IMP.core.Reference.setup_particle(mdl, pi_constr, pi_ref)
        lcp.add(pi_constr)

c = IMP.container.SingletonsConstraint(smp, None, lcp)
mdl.add_score_state(c)
mdl.update()

dof.disable_movers(ps_constr,[IMP.core.BallMover])




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

for csv in ["../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFG.csv","../data/crosslinks/old/Exp2-KFG-081016-2.csv",
            "../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFGE.csv","../data/crosslinks/old/Exp2-KFGE-081016-2.csv",
           "../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFGV.csv","../data/crosslinks/old/Exp3-KFG-081016-3.csv",
            "../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFGVE.csv","../data/crosslinks/old/Exp3-KFGE-081016-3.csv",
            "../data/crosslinks/old/Vgrg-Tle1.csv","../data/crosslinks/old/Vgrg-Tle1-2-1008.csv"]:

#for csv in ["../data/crosslinks/old/Exp1-BP-CXMS-12aout16-KFG.csv","../data/crosslinks/old/Exp2-KFG-081016-2.csv"]:

    cldb_tmp=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    
    cldb_tmp.create_set_from_file(csv)
    
    cldb_tmp.rename_proteins({"TssKHis6":"TssK","StrepTssF":"TssF","TssGFlag":"TssG"})
    cldb_tmp.set_name(csv)
    cldb.append_database(cldb_tmp)
    
    
for xl in cldb:
    xl[cldb.id_score_key]=-math.log(xl[cldb.id_score_key])

cldb.create_new_keyword(cldb.sigma1_key)
cldb.create_new_keyword(cldb.sigma2_key)

cldb.set_value(cldb.sigma1_key,'SIGMA')
cldb.set_value(cldb.sigma2_key,'SIGMA')
cldb.set_value(cldb.sigma1_key,'FLEXIBLE-1',FO(cldb.residue1_links_number_key,operator.gt,5))
cldb.set_value(cldb.sigma2_key,'FLEXIBLE-1',FO(cldb.residue2_links_number_key,operator.gt,5))
cldb.set_value(cldb.sigma1_key,'FLEXIBLE-2',FO(cldb.residue1_links_number_key,operator.gt,30))
cldb.set_value(cldb.sigma2_key,'FLEXIBLE-2',FO(cldb.residue2_links_number_key,operator.gt,30))


#for xl in cldb:
#    print xl[cldb.sigma1_key],xl[cldb.residue1_links_number_key],xl[cldb.sigma2_key],xl[cldb.residue2_links_number_key]


fo=FO("Decision",operator.eq,"Accepted")

fcldb=cldb.filter(fo)

data=IMP.crosslinkms.DataFixedLength("/Bis/home/rpellari/imp-projects/imp-050517/imp/modules/crosslinkms/data/data.21.10000.gaussian.1.rmf")
xlmovers=[]
xls=[]

xl = IMP.crosslinkms.restraint_pmi2.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                CrossLinkDataBase=fcldb,
                                length=21.0,
                                slope=0.01,
                                resolution=1.0,
                                label="XL",
                                data=data,
                                alpha=10,
                                beta=100,
                                isotonic_regression=fcldb.id_score_key,
                                outlier=True) #cldb.id_score_key)
xl.add_to_model()

sigma=xl.sigma_dictionary["SIGMA"][0]
sigma.set_scale(11.0)
xl.set_sigma_is_sampled(True)
sigma=xl.sigma_dictionary["FLEXIBLE-1"][0]
sigma.set_scale(11.0)
sigma=xl.sigma_dictionary["FLEXIBLE-2"][0]
sigma.set_scale(11.0)
xl.set_sigma_is_sampled(True)
psi=xl.psi_dictionary["PSI"][0]
psi.set_scale(0.05)
xl.set_psi_is_sampled(False)
output_objects.append(xl)
xls.append(xl)
xlmovers+=xl.get_movers()




##################  EM Restrains ##############
densities=[]
densities+= IMP.atom.Selection(hier,molecule="TssF",representation_type=IMP.atom.DENSITIES).get_selected_particles()
densities+= IMP.atom.Selection(hier,molecule="TssG",representation_type=IMP.atom.DENSITIES).get_selected_particles()
densities+= IMP.atom.Selection(hier,molecule="TssK",representation_type=IMP.atom.DENSITIES).get_selected_particles()

m=0
for i in densities:
    m+=IMP.atom.Mass(i).get_mass()
print 'moleculer mass',m


gemt = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                            './256.gmm.impgmm',
                                            slope=0.01,
                                            target_radii_scale=3.0,
                                            weight=1.0,
                                            target_mass=448158.5538,
                                            target_is_rigid_body=False)

gemt.add_to_model()
output_objects.append(gemt)
#gemt.write_target_gmm_to_mrc(fileout="output.mrc", voxel_size=5.0)
#exit()
#gemt.center_on_target_density()

#######################################


#IMP.pmi.tools.shuffle_configuration(bs.get_molecules()[0]['TssF'],max_translation=100,avoidcollision_rb=True, avoidcollision_fb=True)
#IMP.pmi.tools.shuffle_configuration(bs.get_molecules()[0]['TssG'],max_translation=100,avoidcollision_rb=True, avoidcollision_fb=True)


cg=IMP.pmi.samplers.ConjugateGradients(mdl,[])
cg.optimize(1000)
cg.optimize(1000)


dof.enable_all_movers()

for xyz in fixedxyz:
    IMP.core.XYZ(xyz).set_coordinates_are_optimized(True)
    
for rb in fixedrbs:
    IMP.core.RigidBody(rb).set_coordinates_are_optimized(True)



init=True
rex_obj=None

import IMP.rmf
import RMF
import math

for k in range(100):

	 rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                        	     root_hier=hier,                          # pass the root hierarchy
                                        	     #crosslink_restraints=xls,                     # will display like XLs
                                        	     monte_carlo_sample_objects=dof.get_movers(),#+xlmovers,# pass MC movers
                                        	     global_output_directory='output_debug_'+str(k),
                                        	     output_objects=output_objects,
                                        	     monte_carlo_steps=10,
                                        	     monte_carlo_temperature=1.0,
                                        	     simulated_annealing=False,
                                        	     simulated_annealing_minimum_temperature=1.0,
                                        	     simulated_annealing_maximum_temperature=2.5,
                                        	     simulated_annealing_minimum_temperature_nframes=100,
                                        	     simulated_annealing_maximum_temperature_nframes=100,
                                        	     replica_exchange_minimum_temperature=1.0,
                                        	     replica_exchange_maximum_temperature=1000.0,
                                        	     number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                        	     number_of_frames=10,
						     replica_exchange_object=rex_obj)                   # increase number of frames to get better results!
	 rex.execute_macro()
         rex_obj=rex.replica_exchange_object
	 

	 gemt.gaussianEM_restraint.debug()
	 indexes=gemt.gaussianEM_restraint.get_indexes()
	 scores=gemt.gaussianEM_restraint.get_log2()

	 scores=[math.log(s) for s in scores]

         if init:


		 map_p={}
		 root_n=IMP.atom.Hierarchy(IMP.Particle(mdl))
		 for n,i in enumerate(indexes):
		     pi=mdl.get_particle(i)
		     coor=IMP.core.XYZ(pi).get_coordinates()
		     p=IMP.Particle(mdl)
		     d=IMP.core.XYZR.setup_particle(p)
		     d.set_coordinates(coor)
		     d.set_radius(5.0)
		     IMP.atom.Mass.setup_particle(p,1.0)
		     col=IMP.display.get_rgb_color((scores[n]-min(scores))/(max(scores)-min(scores)))
		     c=IMP.display.Colored.setup_particle(p,col)
		     h=IMP.atom.Hierarchy.setup_particle(p)
		     p.set_name(str(scores[n]))
		     root_n.add_child(h)
		     map_p[i]=p

		 rh = RMF.create_rmf_file("debug"+str(k)+".rmf")
		 IMP.rmf.add_hierarchies(rh, [root_n])
		 IMP.rmf.save_frame(rh)
		 del rh
                 init=False
		 
	 else:
		for n,i in enumerate(indexes):
		     col=IMP.display.get_rgb_color((scores[n]-min(scores))/(max(scores)-min(scores)))
		     c=IMP.display.Colored(map_p[i]).set_color(col)
		     h=IMP.atom.Hierarchy.setup_particle(p)
		rh = RMF.create_rmf_file("debug"+str(k)+".rmf")
		IMP.rmf.add_hierarchies(rh, [root_n])
		IMP.rmf.save_frame(rh)
		del rh


exit()




print(sf.evaluate(False))
cg.optimize(1000)
print(sf.evaluate(False))
cg.optimize(1000)
print(sf.evaluate(False))





output_dir='./'
output_prefix=''
## Run replica exchange Monte Carlo sampling
rex_obj=None
weights = 1.0 * np.array([1, 1]*100)
i=0;
for w in weights:
    #gemt.set_weight(w)
    print rex_obj
    num_frames=10
    if w==1:
        num_frames=10
    if rex_obj is None:
        rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                            root_hier=hier,                          # pass the root hierarchy
                                            #crosslink_restraints=xls,                     # will display like XLs
                                            monte_carlo_sample_objects=dof.get_movers(),#+xlmovers,# pass MC movers
                                            global_output_directory='%s/%soutput_%d_%g/'%(output_dir,output_prefix, i,w),
                                            output_objects=output_objects,
                                            monte_carlo_steps=10,
                                            monte_carlo_temperature=1.0,
                                            simulated_annealing=False,
                                            simulated_annealing_minimum_temperature=1.0,
                                            simulated_annealing_maximum_temperature=2.5,
                                            simulated_annealing_minimum_temperature_nframes=100,
                                            simulated_annealing_maximum_temperature_nframes=100,
                                            replica_exchange_minimum_temperature=1.0,
                                            replica_exchange_maximum_temperature=1000.0,
                                            number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                            number_of_frames=num_frames)                   # increase number of frames to get better results!
        rex.execute_macro()
        rex_obj=rex.replica_exchange_object
        print dof.rigid_bodies
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
