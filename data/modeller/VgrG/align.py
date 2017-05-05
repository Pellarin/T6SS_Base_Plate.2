from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='../../pdb/4mtk_VgrG.pdb', model_segment=('FIRST:A','LAST:A'))


print mdl


aln.append_model(mdl, align_codes='4mtkA', atom_files='../../pdb/4mtk_VgrG.pdb')
aln.append(file='vgrg.ali', align_codes='VGRG')
aln.align2d(local_alignment=True)
aln.write(file='4mtkA_VgrG.ali', alignment_format='PIR')
aln.write(file='4mtkA_VgrG.pap', alignment_format='PAP')
