import IMP
import IMP.core
import math
import IMP.pmi
import IMP.pmi.tools
import RMF
import IMP.rmf

def constrain_symmetry(mdl,references,
                       clones,
                       transform):

    # get all RBs and particles

    ref_rbs,ref_beads = IMP.pmi.tools.get_rbs_and_beads(references)
    clones_rbs,clones_beads = IMP.pmi.tools.get_rbs_and_beads(clones)

    # dumb check for matching numbers of particles
    if len(ref_rbs)!=len(clones_rbs) or len(ref_beads)!=len(clones_beads):
        raise Exception("ERROR: Your references don't match your clones")

    # symmetry RBs
    for ref,clone in zip(ref_rbs+ref_beads,clones_rbs+clones_beads):
        IMP.core.Reference.setup_particle(clone,ref)
    sm = IMP.core.TransformationSymmetry(transform.get_inverse())
    lsc = IMP.container.ListSingletonContainer(
        mdl,[p.get_particle_index() for p in clones_rbs+clones_beads])
    c = IMP.container.SingletonsConstraint(sm, None, lsc)
    mdl.add_score_state(c)
    print('Created symmetry restraint for',len(ref_rbs),'rigid bodies and',
          len(ref_beads),'flexible beads')
    mdl.update()
    
    
m=IMP.Model()
p=IMP.Particle(m)
h=IMP.atom.Hierarchy.setup_particle(p)

clones=[]
nclones=6
rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0)

for n in range(nclones):
    p0=IMP.Particle(m)
    p1=IMP.Particle(m)
    p2=IMP.Particle(m)
    p3=IMP.Particle(m)

    d0=IMP.core.XYZR.setup_particle(p0)
    d1=IMP.core.XYZR.setup_particle(p1)
    d2=IMP.core.XYZR.setup_particle(p2)
    d3=IMP.core.XYZR.setup_particle(p3)

    d0.set_coordinates((0,0,0))
    d1.set_coordinates((1,0,0))
    d2.set_coordinates((0,1,0))
    d3.set_coordinates((0,0,1))

    d0.set_radius(0.5)
    d1.set_radius(0.5)
    d2.set_radius(0.5)
    d3.set_radius(0.5)
    
    h0=IMP.atom.Hierarchy.setup_particle(p0)
    h1=IMP.atom.Hierarchy.setup_particle(p1)
    h2=IMP.atom.Hierarchy.setup_particle(p2)
    h3=IMP.atom.Hierarchy.setup_particle(p3)

    h0=IMP.atom.Mass.setup_particle(p0,1)
    h1=IMP.atom.Mass.setup_particle(p1,1)
    h2=IMP.atom.Mass.setup_particle(p2,1)
    h3=IMP.atom.Mass.setup_particle(p3,1)
    
    h.add_child(h0)
    h.add_child(h1)
    h.add_child(h2)
    h.add_child(h3)
    
    IMP.atom.create_rigid_body([p0,p1,p2,p3])
    clones.append([p0,p1,p2,p3])
    
    
for n,c in enumerate(clones[1:]):
    rotation_angle = 2.0 * math.pi / float(nclones) * float(n+1)
    rotation3D = IMP.algebra.get_rotation_about_axis(rotational_axis, rotation_angle)
    transformation3D =IMP.algebra.get_rotation_about_point(IMP.algebra.Vector3D(3, 3,   0.00000000), rotation3D)

    constrain_symmetry(m,
                       clones[0],
                       c,
                       transformation3D)

rh = RMF.create_rmf_file("test_symmetry.rmf")
IMP.rmf.add_hierarchies(rh, [h])
IMP.rmf.save_frame(rh)
del rh
