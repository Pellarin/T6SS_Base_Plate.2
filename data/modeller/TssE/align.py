from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='../../pdb/5iw9_TssE.pdb', model_segment=('FIRST:A','LAST:A'))


print mdl


aln.append_model(mdl, align_codes='5iw9A', atom_files='../../pdb/5iw9_TssE.pdb')
aln.append(file='tsse.ali', align_codes='TSSE')
aln.align2d(local_alignment=True)
aln.write(file='5iw9A_TssE.ali', alignment_format='PIR')
aln.write(file='5iw9A_TssE.pap', alignment_format='PAP')
