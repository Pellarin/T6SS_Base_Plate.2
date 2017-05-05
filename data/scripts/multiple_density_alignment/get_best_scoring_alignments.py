import IMP
import IMP.pmi
import IMP.pmi.output
import IMP.pmi.io
import glob
import RMF
import IMP.rmf



best_models = IMP.pmi.io.get_best_models(glob.glob("output_shuffle_rex_xl_test_2/stat.*.out"),
                                                     "Total_Score",
                                                     ["GaussianEMRestraint_sigma_../TssK.50.txt","GaussianEMRestraint_sigma_../TssKFG.256.txt","GaussianEMRestraint_sigma_../TssKFG.256.txt../TssK.50.txt"],
                                                     "rmf_file",
                                                     "rmf_frame_index",
                                                     0,
                                                     1)

rmf_file_list=best_models[0]
rmf_file_frame_list=best_models[1]
score_list=[float(a) for a in best_models[2]]
feature_keyword_list_dict=best_models[3]


#IMP.isd.gmm_tools.write_gmm_to_map(pss,"clusters/"+str(i)+".mrc",5.0,bounding_box=None,origin=None, fast=True)

minindex=score_list.index(min(score_list))

print "best rmf %s best frame %s best score %s"  % (rmf_file_list[minindex],rmf_file_frame_list[minindex],str(score_list[minindex]))
for k in feature_keyword_list_dict:
    print k,feature_keyword_list_dict[k][minindex]


model=IMP.Model()
rh = RMF.open_rmf_file_read_only(rmf_file_list[minindex])
prots = IMP.rmf.create_hierarchies(rh, model)
IMP.rmf.load_frame(rh, RMF.FrameID(rmf_file_frame_list[minindex]))

for n,p in enumerate(prots[0].get_children()):
    ps=IMP.atom.get_leaves(p)
    pss=[IMP.core.Gaussian(p) for p in ps]
    IMP.isd.gmm_tools.write_gmm_to_map(pss,str(n)+".mrc",5.0,bounding_box=None,origin=None, fast=True)
