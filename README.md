# T6SS_Base_Plate

Information

For VgrG the closest homolog and more complete structure will be pdb_id 4MTK

For TssE:21-128 | pdb 2ia7_A, is OK. But there is a very recent pdb modeled in Taylor… Leiman PG, Nature 2016, 4hrz

For TssA, I solved two X-ray structure: TssANt and TssACt

TssA:221-377 | pdb 4yo3_A; Identities=100%  Similarity=1.519

TssA:406-529 | pdb 4yo5_A; Identities=100%  Similarity=1.519

TssK: only SAXS and EM-NS low resolution models

TssF and TssG are distantly related to two T4 bacteriophages proteins (small domains conserved; cf Brunet YR… Cascales, PloS Genet 2015)

To do:

  1) An homology model of VgrG
  
  2) Homology models for TssF and TssG, if possible
  
  3) Homology model of TssE, based on Taylor… Leiman PG, Nature 2016

1) Composition/Stoichiometry. Discrepnacy between biochemal staining and native MS. For instance: KFG 2F+5K+1G by stainng, 6K+2F+1G by mass spec. K is trimer, alone, and this makes sense with the EM. TssE is probably stochiometry of 1. VgrG is a Trimer, we know that.
2) EM of subcomplexes ( no symmetry ). 28-30 angstrom. Saxs envelope o
f K and em model of K. Conflict between the density of KFG and the predicted stoichiomtry.
4) pull K and get all proteins.
3) Symmetry. No symmetry.
5) Ion collision NatvMS remove first TssE, (TssE, KFG) 



Urgent: KFG and KFGE, we have EM negative-staining. Redo the XL data, new dataset of XL to increase the confidence. Directed mutagenesis to see the effect of mutations on cross-linked spots, you will see on the offect in vivo (assembly and activity). Full complex: KFGEV,we don't have the EM model.








Dear Riccardo,

Please find a summary of published (Flaugnatti et al. Mol Microbiol 2016, attached to this email) and unpublished results that could be useful, in addition to the EM data. 

-Tle1 is an anti-bacterial effector with phospholipase A1 activity (the catalytic S297 is required for PLA1 activity. The 2 other putative residues of the catalytic triad: Asp245, His310). Note that Nicolas has evidence that Tle1 is not active (no PLA1 activity) when in complex with VgrG.

-Tle1 is encoded downstream the vgrG gene in the T6SS-1 gene cluster. VgrG and Tle1 Stop and Start codons overlap.

-We showed the interaction of Tle1 with VgrG by co-immunoprecipitation (coIP) and bacterial two-hybrid (BACTH) and defined regions in VgrG necessary and sufficient for this interaction. 

-We divided VgrG in three domains: VgrG-DUF2345-TTR. VgrG1 carries an additional C-terminal domain (CTD) separated from the gp27-gp5 common core by a predicted coiled-coil region. VgrG CTD amino-acid 611-766 region is predicted (by HHPred or Phyre2) to be constituted of a regular repetition of small short-strands that are reminiscent to the C-terminal domain of gp5 and likely extends the VgrG spike. This additional b-prism domain is followed by a 62-amino-acid region (residues 780 to 841) predicted to fold as a transthyretin (TTR)-like domain. 

-Regions in VgrG: (Flaugnatti et al. , Fig7)
From coIP and BACTH analyses using VgrG truncated derivatives, we concluded that Tle1 interacts with VgrG C-terminal motif (aa 771-841), and maybe a second motif is present in the C-ter region of VgrG (aa 616-841).

-Regions in Tle1 involved in Tle1/VgrG interaction: (unpublished)

Upon alignment of Tle1 orthologs, we noticed that they differ in their N-terminal part: Tle1 orthologs possess either a N-terminal extension of about 30 amino-acids, or a region encoding a Hcp protein (or a PAAR containing protein). As Hcp and PAAR were shown or proposed to be used for toxin transport, we thought that this N-terminal 30 a.a extension of Tle1could be required for Tle1 export (and the other orthologs could have another mode of transport).
We thus deleted the first 26 aa of Tle1 in coIP experiments, and Tle1(delta1-26) interaction with VgrG was strongly affected by the deletion (but not completely abolished).
(Note that the resulting truncated derivative seems to be still active as it retains its toxicity when produce in the periplasm of E. coli). 
-In vivo, all the mutants we identified (in vgrG or tle1) that abolish the interaction are unable to perform T6SS antibacterial killing.   
In conclusion, in our hands, the C-terminal TTR domain of VgrG and the N-terminal domain of Tle1 are required for the interaction, but may be not the only determinants of the interaction (as deletion of these two domains strongly affect but does not completely abolish the interaction by coIP).
These results are in good agreement with the cross-links found by Martial. We now have to test the other regions identified. This is where we are at the moment!
Happy Christmas et  bonnes fêtes à tous!
Cheers,
Laure
