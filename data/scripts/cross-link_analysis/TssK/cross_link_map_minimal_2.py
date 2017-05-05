import IMP
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.tools
import IMP.pmi.io.xltable
import IMP.pmi.restraints.crosslinking_new
from IMP.pmi.io.crosslink import FilterOperator as FO
import operator
import math

fastadirectory="./"
xlmsdirectory="./"

rplp=IMP.pmi.io.crosslink.ResiduePairListParser("MSSTUDIO")
cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter(rplp)
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_site_pairs_key("Selected Sites")
cldbkc.set_id_score_key("E")
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.set_name("master")



for csv in ["../../../crosslinks/Exp1-BP-CXMS-12aout16-KFG.csv","../../../crosslinks/Exp2-KFG-081016-2.csv",
            "../../../crosslinks/Exp1-BP-CXMS-12aout16-KFGE.csv","../../../crosslinks/Exp2-KFGE-081016-2.csv",
            "../../../crosslinks/Exp1-BP-CXMS-12aout16-KFGV.csv","../../../crosslinks/Exp3-KFG-081016-3.csv",
            "../../../crosslinks/Exp1-BP-CXMS-12aout16-KFGVE.csv","../../../crosslinks/Exp3-KFGE-081016-3.csv",
            "../../../crosslinks/Vgrg-Tle1.csv","../../../crosslinks/Vgrg-Tle1-2-1008.csv"]:

    cldb_tmp=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    
    cldb_tmp.create_set_from_file(csv)
    
    cldb_tmp.rename_proteins({"TssKHis6":"TssK"})
    cldb_tmp.set_name(csv)
    cldb.append_database(cldb_tmp)


fo=FO("Decision",operator.eq,"Accepted")

fcldb=cldb.filter(fo)


prots =    ["TssK"]

protx=prots 
proty=prots

suffix="TssK"

conversion_fastanames_name = {"TssKHis6":"TssK"}
                           
inv_conversion_fastanames_name = {v: k for k, v in conversion_fastanames_name.items()}

sequences={}
for i  in prots:
    sequences[i]=("../../../Basteplate.crosslink.fasta",  inv_conversion_fastanames_name[i], i)

xlt1=IMP.pmi.io.xltable.XLTable(20)
for p in list(IMP.pmi.tools.OrderedSet(protx+proty)):
    t=sequences[p]
    xlt1.load_sequence_from_fasta_file(t[0],  t[1], t[2]) 

xlt1.load_rmf_coordinates("output_"+suffix+"/rmfs/0.rmf3",99,list(IMP.pmi.tools.OrderedSet(protx+proty)),nomap=True)

#dictionary_of_gaps={}

xlt1.load_crosslinks(fcldb)
xlt1.compute_distances()
fcldb.save_csv("xl_distances_"+suffix+".csv")
fcldb.plot("output.pdf",type="scatter",xkey="MinAmbiguousDistance",ykey=cldbkc.redundancy_key,colorkey=cldbkc.id_score_key)
xlt1.save_rmf_snapshot("bestsnapshot_"+suffix+".rmf",color_id=cldbkc.redundancy_key)

