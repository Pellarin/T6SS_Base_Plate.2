import IMP
import IMP.container
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import RMF
import math


mdl=IMP.Model()
root=IMP.atom.Hierarchy(IMP.Particle(mdl))

rbs_ps={}
ps_mv={}
mvs=[]
pss=[]
for mol in range(10):
    mol=IMP.atom.Hierarchy(IMP.Particle(mdl))
    root.add_child(mol)
    ps=[]
    for c in [(0,0,0),(5,0,0),(0,5,0)]:
        p=IMP.Particle(mdl)
        dr=IMP.core.XYZR.setup_particle(p)
        dr.set_coordinates(c)
        dr.set_radius(1.0)
        IMP.atom.Mass.setup_particle(p,1.0)
        h=IMP.atom.Hierarchy(p)
        mol.add_child(h)
        ps.append(p)
    rb = IMP.atom.create_rigid_body(IMP.atom.get_leaves(mol))
    rb_mover = IMP.core.RigidBodyMover(rb.get_model(),rb,1,0.5)
    mvs.append(rb_mover)
    rbs_ps[rb]=ps
    pss.append(ps)
    for p in ps:
        if p in ps_mv:
            ps_mv[p].append(rb_mover)
        else:
            ps_mv[p]=[rb_mover]



#setup 2 floppy bodies

floatkeys = [IMP.FloatKey(4), IMP.FloatKey(5), IMP.FloatKey(6)]
for rb in rbs_ps:
    for p in rbs_ps[rb][1:2]:
        rb.set_is_rigid_member(p.get_index(),False)
        for fk in floatkeys:
            p.set_is_optimized(fk,True)
        fbmv=IMP.core.BallMover(p.get_model(), p,
                                            IMP.FloatKeys(floatkeys),
                                            1.0)
        mvs.append(fbmv)

        if p in ps_mv:
            ps_mv[p].append(fbmv)
        else:
            ps_mv[p]=[fbmv]


# get all RBs and particles
def add_symmetry(clones,refs,transform):
    href    = IMP.pmi.tools.input_adaptor(refs,'all',flatten=True)
    hclones = IMP.pmi.tools.input_adaptor(clones,'all',flatten=True)

    ref_rbs,ref_beads = IMP.pmi.tools.get_rbs_and_beads(href)
    clones_rbs,clones_beads = IMP.pmi.tools.get_rbs_and_beads(hclones)

    for ref,clone in zip(ref_rbs,clones_rbs):
        IMP.core.Reference.setup_particle(clone,ref)
    
    for ref,clone in zip(ref_beads,clones_beads):
        IMP.core.Reference.setup_particle(clone,ref)
    
    sm = IMP.core.TransformationSymmetry(transform)
    lsc = IMP.container.ListSingletonContainer(mdl,[p.get_particle().get_index() for p in clones_rbs+clones_beads])
    c = IMP.container.SingletonsConstraint(sm, None, lsc)

    mdl.add_score_state(c)
        
rotational_axis=IMP.algebra.Vector3D(-0.98477941,0.15168411,-0.08486137)
rotational_point=IMP.algebra.Vector3D(0.00000000,0,0)

disabled_movers=set()
for n,c in enumerate(pss[1:]):
    print(n,c)
    rotation_angle = 2.0 * math.pi / len(pss) * float(n+1)
    rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
    transformation3D =IMP.algebra.get_rotation_about_point(rotational_point, rotation3D)
    add_symmetry(c,pss[0],transformation3D)
    for p in c:
        for mv in ps_mv[p]:
            disabled_movers.add(mv)
    
mdl.update()

enabled_movers=[]
for mv in mvs:
    if mv in disabled_movers:
        continue
    else:
        enabled_movers.append(mv)


print enabled_movers

ssps = IMP.core.SoftSpherePairScore(1.0)
lsa = IMP.container.ListSingletonContainer(mdl)
lsa.add(IMP.get_indexes(IMP.atom.get_leaves(root)))
rbcpf = IMP.core.RigidClosePairsFinder()
cpc = IMP.container.ClosePairContainer(lsa, 0.0, rbcpf, 10.0)
evr = IMP.container.PairsRestraint(ssps, cpc)

sf = IMP.core.RestraintsScoringFunction([evr])
mc = IMP.core.MonteCarlo(mdl)
mc.set_scoring_function(sf)

sm = IMP.core.SerialMover(enabled_movers)
mc.add_mover(sm)
mc.set_return_best(False)
mc.set_kt(1.0)

rh = RMF.create_rmf_file("out.rmf")
IMP.rmf.add_hierarchies(rh, [root])

for nloop in range(10000):
    mc.optimize(10)
    IMP.rmf.save_frame(rh)
    rh.flush()

del rh