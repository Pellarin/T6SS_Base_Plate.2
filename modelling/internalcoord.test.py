import IMP
import IMP.container
import IMP.atom
import IMP.rmf
import IMP.pmi
import RMF


mdl=IMP.Model()
root=IMP.atom.Hierarchy(IMP.Particle(mdl))
mol1=IMP.atom.Hierarchy(IMP.Particle(mdl))
mol2=IMP.atom.Hierarchy(IMP.Particle(mdl))
root.add_child(mol1)
root.add_child(mol2)

rbs_ps={}
mvs=[]
for mol in [mol1,mol2]:
    ps=[]
    for c in [(0,0,0),(1,0,0),(0,1,0)]:
        p=IMP.Particle(mdl)
        dr=IMP.core.XYZR.setup_particle(p)
        dr.set_coordinates(c)
        dr.set_radius(1.0)
        IMP.atom.Mass.setup_particle(p,1.0)
        h=IMP.atom.Hierarchy(p)
        mol.add_child(h)
        ps.append(p)
    rb = IMP.atom.create_rigid_body(IMP.atom.get_leaves(mol))
    rb_mover = IMP.core.RigidBodyMover(rb.get_model(), rb, 1,
                                           0.5)
    mvs.append(rb_mover)
    rbs_ps[rb]=ps
    
#setup 2 floppy bodies

floatkeys = [IMP.FloatKey(4), IMP.FloatKey(5), IMP.FloatKey(6)]
for rb in rbs_ps:
    for p in rbs_ps[rb][0:2]:
        rb.set_is_rigid_member(p.get_index(),False)
        for fk in floatkeys:
            p.set_is_optimized(fk,True)
        fbmv=IMP.core.BallMover(p.get_model(), p,
                                            IMP.FloatKeys(floatkeys),
                                            1.0)
        mvs.append(fbmv)





class _InternalCoordinatesConstraint(IMP.SingletonModifier):

    """A class that updates internal coordinates between rigid bodies
    """

    def __init__(self,m):
        IMP.SingletonModifier.__init__(
            self, "InternalCoordinatesModifier%1%")
        
            
    def apply_index(self, m, pi):
        ref_p=IMP.core.Reference(m,pi).get_reference_particle()
        ref_index=ref_p.get_index()
        nrm1=IMP.core.NonRigidMember(m,ref_index)
        nrm2=IMP.core.NonRigidMember(m,pi)
        c1=nrm1.get_internal_coordinates()
        nrm2.set_internal_coordinates(c1)

    def do_get_inputs(self, m, pis):
        return IMP.get_particles(m, pis)

    def do_get_outputs(self, m, pis):
        return self.do_get_inputs(m, pis)



ps_constr=[]
lcp = IMP.container.ListSingletonContainer(mdl)
smp = IMP.pmi.InternalCoordinatesConstraint() 
#smp = _InternalCoordinatesConstraint(mdl)
for (pref,pconstr) in zip(IMP.atom.get_leaves(mol1),IMP.atom.get_leaves(mol2)):
    pi_ref=pref.get_particle_index()
    pi_constr=pconstr.get_particle_index()
    ps_constr.append(pconstr.get_particle())
    
    if IMP.core.NonRigidMember.get_is_setup(mdl,pi_ref) and IMP.core.NonRigidMember.get_is_setup(mdl,pi_constr):
        
        IMP.core.Reference.setup_particle(mdl, pi_constr, pi_ref)
        lcp.add(pi_constr)

c = IMP.container.SingletonsConstraint(smp, None, lcp)
mdl.add_score_state(c)
mdl.update()



ssps = IMP.core.SoftSpherePairScore(1.0)
lsa = IMP.container.ListSingletonContainer(mdl)
lsa.add(IMP.get_indexes(IMP.atom.get_leaves(root)))
rbcpf = IMP.core.RigidClosePairsFinder()
cpc = IMP.container.ClosePairContainer(lsa, 0.0, rbcpf, 10.0)
evr = IMP.container.PairsRestraint(ssps, cpc)

sf = IMP.core.RestraintsScoringFunction([evr])
mc = IMP.core.MonteCarlo(mdl)
mc.set_scoring_function(sf)

sm = IMP.core.SerialMover(mvs)
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