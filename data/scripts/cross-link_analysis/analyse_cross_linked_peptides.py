import IMP
import IMP.pmi.io
import IMP.pmi.io.crosslink
from IMP.pmi.io.crosslink import FilterOperator as FO
import operator



histo_dict={}
histo_dict_inter={}
histo_dict_pairwise={}
peptide_dict_pairwise={}

seqs=IMP.pmi.topology.Sequences("../../Basteplate.crosslink.fasta")

for prot in seqs:
    histo_dict[prot]=[0]*len(seqs[prot])
    histo_dict_inter[prot]=[0]*len(seqs[prot])
    peptide_dict_pairwise[prot]={}
    for prot1 in seqs:
        histo_dict_pairwise[(prot,prot1)]=[0]*len(seqs[prot])

#for csv in ["../../crosslinks/Exp1-BP-CXMS-12aout16-KFG.csv","../../crosslinks/Exp2-KFG-081016-2.csv",
#            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGE.csv","../../crosslinks/Exp2-KFGE-081016-2.csv",
#            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGV.csv","../../crosslinks/Exp3-KFG-081016-3.csv",
#            "../../crosslinks/Exp1-BP-CXMS-12aout16-KFGVE.csv","../../crosslinks/Exp3-KFGE-081016-3.csv"]:

for csv in [ "../../crosslinks/old/Exp3-KFGE-081016-3.csv"]:


    rplp=IMP.pmi.io.crosslink.ResiduePairListParser("MSSTUDIO")
    cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter(rplp)
    cldbkc.set_protein1_key("Protein 1")
    cldbkc.set_protein2_key("Protein 2")
    cldbkc.set_site_pairs_key("Selected Sites")
    cldbkc.set_id_score_key("E")
    cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb.create_set_from_file(csv)
    nxl=0
    for xl in cldb:
        try:
           peptide1start=int(xl["Start 1"])
           peptide1end=int(xl["Stop 1"])
           peptide2start=int(xl["Start 2"])
           peptide2end=int(xl[" Stop 2"])
           peptide1=xl["Peptide 1"]
           peptide2=xl["Peptide 2"]
           decision=xl["Decision"]
           type=xl["Type"]
           intrainter=xl["Intra/Inter"]
        except:
           continue 
        if not decision == "Accepted":
           continue 
        if not type == "Cross-linked":
           continue
        if not intrainter == "Inter":
           inter=False
        else:
           inter=True 
        
        ambiguity=xl[cldb.ambiguity_key]
        protein1=xl[cldb.protein1_key]
        protein2=xl[cldb.protein2_key]
        residue1=xl[cldb.residue1_key]
        residue2=xl[cldb.residue2_key]
        
        if protein1 == "StrepTssF":
            if residue1 <= 9 and residue1 >= 2:
                continue

        if protein2 == "StrepTssF":
            if residue2 <= 9 and residue2 >= 2:
                continue

        if protein1 == "TssKHis6":
            if residue1 <= 451 and residue1 >= 446:
                continue

        if protein2 == "TssKHis6":
            if residue2 <= 451 and residue2 >= 446:
                continue


        if protein1 == "TssGFlag":
            if residue1 <= 374 and residue1 >= 367:
                continue

        if protein2 == "TssGFlag":
            if residue2 <= 374 and residue2 >= 367:
                continue

        if protein1 == "VgrGHA":
            if residue1 <= 850 and residue1 >= 843:
                continue

        if protein2 == "VgrGHA":
            if residue2 <= 850 and residue2 >= 843:
                continue

        if protein1 == "HATssE":
            if residue1 <= 10 and residue1 >= 2:
                continue

        if protein2 == "HATssE":
            if residue2 <= 10 and residue2 >= 2:
                continue
        
        if protein1 and protein2:
           nxl+=1.0/ambiguity
           list1=[]
           for i in range(1,peptide1end+1):
               if i<peptide1start:
                  list1.append(0)
               else:
                  list1.append(1.0/ambiguity)
           
           list2=[]
           for i in range(1,peptide2end+1):
               if i<peptide2start:
                  list2.append(0)
               else:
                  list2.append(1.0/ambiguity)
           
           nmax=len(histo_dict_pairwise[(protein1,protein2)])
           for n,i in enumerate(list1):
               if n<nmax: 
                  histo_dict_pairwise[(protein1,protein2)][n]+=i
               else:
                  histo_dict_pairwise[(protein1,protein2)].append(i)  
                                 
           nmax=len(histo_dict_pairwise[(protein2,protein1)])
           for n,i in enumerate(list2):
               if n<nmax: 
                  histo_dict_pairwise[(protein2,protein1)][n]+=i
               else:
                  histo_dict_pairwise[(protein2,protein1)].append(i) 

           crosslinked_peptides12=(protein1,peptide1start,peptide1,peptide1end,protein2,peptide2start,peptide2,peptide2end)
           crosslinked_peptides21=(protein2,peptide2start,peptide2,peptide2end,protein1,peptide1start,peptide1,peptide1end)
           
           if not crosslinked_peptides12 in peptide_dict_pairwise[protein1]:
              peptide_dict_pairwise[protein1][crosslinked_peptides12]=1.0/ambiguity
           else:
              peptide_dict_pairwise[protein1][crosslinked_peptides12]+=1.0/ambiguity
              
           if not crosslinked_peptides21 in peptide_dict_pairwise[protein2]:
              peptide_dict_pairwise[protein2][crosslinked_peptides21]=1.0/ambiguity
           else:
              peptide_dict_pairwise[protein2][crosslinked_peptides21]+=1.0/ambiguity           

           
           for protein,list in [(protein1,list1),(protein2,list2)]:

               nmax=len(histo_dict[protein])
               for n,i in enumerate(list):
                   if n<nmax: 
                      histo_dict[protein][n]+=i
                   else:
                      histo_dict[protein].append(i)                   
               
               if inter:
                   nmax=len(histo_dict_inter[protein])
                   for n,i in enumerate(list):
                       if n<nmax: 
                          histo_dict_inter[protein][n]+=i
                       else:
                          histo_dict_inter[protein].append(i)                    
                

        

       
#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

for prot in peptide_dict_pairwise:
    for pept in peptide_dict_pairwise[prot]:
        print peptide_dict_pairwise[prot][pept],pept[0]+" "+str(pept[1])+" "+pept[2]+" "+str(pept[3])+" "+pept[4]+" "+str(pept[5])+" "+pept[6]+" "+str(pept[7])
        

pnames=['StrepTssF', 'HATssE', 'VgrGHA', 'TssKHis6', 'TssGFlag']

handles=[]
for protein1 in pnames:
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    h,=plt.plot(range(1,len(histo_dict[protein1])+1),histo_dict[protein1],":b",label="total",linewidth=2.0)
    handles.append(h)
    for protein2 in pnames:
        h,=plt.plot(range(1,len(histo_dict_pairwise[(protein1,protein2)])+1),histo_dict_pairwise[(protein1,protein2)],label=protein2)
        handles.append(h)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)
    ax.set_xlabel(protein1)
    plt.savefig(protein1+"KFGE.3.pdf", format='pdf')
    plt.show()
               