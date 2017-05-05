import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking_new
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.io

import IMP.bayesianem
import IMP.bayesianem.restraint

import os
import operator
import math

class _NuisanceSingletonModifier(IMP.SingletonModifier):

    """A class that updates the membrane particles
    """

    def __init__(self,sigma1,sigma2,sigma3):
        IMP.SingletonModifier.__init__(
            self, "NuisanceSingletonModifier%1%")
        self.sigma1=sigma1
        self.sigma2=sigma2
        self.sigma3=sigma3

    def apply_index(self, m, pi):
        sigma2_value = self.sigma2.get_nuisance()
        sigma3_value = self.sigma3.get_nuisance()
        d = IMP.isd.Nuisance(m, pi)
        d.set_nuisance(sigma3_value/sigma2_value)

    def do_get_inputs(self, m, pis):
        return IMP.get_particles(m, pis)

    def do_get_outputs(self, m, pis):
        return self.do_get_inputs(m, pis)


def create_gem_ref(fix_ps,filename):
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(fix_ps,filename,
                                                slope=0.000001,
                                                target_radii_scale=3.0,
                                                target_is_rigid_body=True)
    gem.set_label(filename)
    gem.add_target_density_to_hierarchy(root_hier)
    gem.add_to_model()
    
    rb=gem.rb
    rbxyz = (rb.get_x(), rb.get_y(), rb.get_z())

    transformation = IMP.algebra.get_random_local_transformation(
        rbxyz,
        100,
        math.pi)

    IMP.core.transform(rb, transformation)
    
    return gem

def create_gem_cross(gem1,gem2):

    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(gem1.target_ps,target_ps=gem2.target_ps,
                                                slope=0.000001,
                                                target_radii_scale=3.0,
                                                target_is_rigid_body=False)
    gem.set_label(gem1.label+gem2.label)

    smp = _NuisanceSingletonModifier(gem.sigmaglobal,gem1.sigmaglobal,gem2.sigmaglobal)
    lcp = IMP.container.ListSingletonContainer(m)
    lcp.add(gem.sigmaglobal.get_particle().get_index())
    c = IMP.container.SingletonsConstraint(smp, None, lcp)
    m.add_score_state(c)
    gem.sigmaissampled=False

    gem.add_to_model()
    return gem


# setting up parameters

rbmaxtrans = 4.00
fbmaxtrans = 5.00
rbmaxrot=0.05
outputobjects = []
sampleobjects = []

# setting up topology

m = IMP.Model()
s = IMP.pmi.topology.System(m)
root_hier = s.build()

em_directory="../../em/"

#em
fixed_ps = []
IMP.isd.gmm_tools.decorate_gmm_from_text(
    '../TssKFGE.256.txt',
    fixed_ps,
    m,
    radius_scale=3.0,
    mass_scale=1.0)
    
f=IMP.atom.Fragment().setup_particle(IMP.Particle(m))
f.set_name("GaussianEMRestraint_density_Target")
for p in fixed_ps:
  f.add_child(p)
root_hier.add_child(f)

gems=[]
for gmm in ["../TssKFG.256.txt","../TssK.50.txt"]:
    gem=create_gem_ref(fixed_ps,gmm)
    outputobjects.append(gem)
    sampleobjects.append(gem)
    gems.append(gem)
    print(gem.get_particles_to_sample())
    print(gem.get_output())

gem=create_gem_cross(gems[0],gems[1])
outputobjects.append(gem)
sampleobjects.append(gem)


mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hier,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=None,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=20.0,
                                    simulated_annealing_minimum_temperature_nframes=100,
                                    simulated_annealing_maximum_temperature_nframes=100,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=100.0,
                                    number_of_best_scoring_models=0,
                                    monte_carlo_steps=10,
                                    number_of_frames=100000,
                                    write_initial_rmf=True,
                                    initial_rmf_name_suffix="initial",
                                    stat_file_name_suffix="stat",
                                    best_pdb_name_suffix="model",
                                    do_clean_first=True,
                                    do_create_directories=True,
                                    global_output_directory="output_shuffle_rex_xl_test_2",
                                    rmf_dir="rmfs/",
                                    best_pdb_dir="pdbs/",
                                    replica_stat_file_suffix="stat_replica")
mc1.execute_macro()

