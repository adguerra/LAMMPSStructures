import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation


sim = Simulation("run", 3, 0.1, 0.03, 0.1)

beam_length = 0.095
beam_stiffness = 10 ** 5
beam_thickness = 0.003

grain_stiffness = 10 ** 5
grain_diameter = 0.015

timestep = 0.0000001
density = 0.01
viscosity = 0.000001

# We need to manually edit the timestep because we are about to both connect_particles_to_elastic (which executes a run 0), 
# AND move (which requres that the timestep be unchanged between runs)
sim.manually_edit_timestep(timestep)

beam_type, bond_type, angle_type = sim.add_beam(50, np.array([0, 0, -beam_length/2]),np.array([0, 0, beam_length/2]), beam_thickness, beam_stiffness, density)
grain_type = sim.add_grains(np.array([[-grain_diameter/2 - beam_thickness/2, 0, grain_diameter/2], [-grain_diameter/2 - beam_thickness/2, 0, -grain_diameter/2]]), grain_diameter, density)
grain_type2 = sim.add_grains(np.array([[grain_diameter/2 + beam_thickness/2, 0, grain_diameter/2], [grain_diameter/2 + beam_thickness/2, 0, -grain_diameter/2]]), grain_diameter, density)

sim.connect_particles_to_elastic(beam_type, 1, type = grain_type, n_bonds_per_grain=2)
sim.connect_particles_to_elastic(beam_type, 1, type = grain_type2, n_bonds_per_grain=2)

sim.move(particles = 1, xvel = 0, yvel = 0, zvel = 0.01)
sim.move(particles = 50, xvel = 0, yvel = 0, zvel = -0.01)
sim.perturb(particles = 25, xdir = 1)

sim.turn_on_granular_potential(grain_type, grain_type, grain_stiffness)
sim.turn_on_granular_potential(grain_type2, grain_type2, grain_stiffness)

sim.add_viscosity(viscosity)
sim.design_dump_files(0.05)
sim.run_simulation(1)




