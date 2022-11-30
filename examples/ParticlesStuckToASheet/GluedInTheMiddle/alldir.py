import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":
    simname = str(sys.argv[1])
    

    # Start a simulation with the name simname
    sim = Simulation(simname, 3, 0.5, 0.5, 0.2)

    grain_diameter = 0.03
    sheet_radius = 0.159
    grain_packing_fraction = 1
    number_of_grains = int(np.round((sheet_radius ** 2) * grain_packing_fraction / ((grain_diameter / 2) ** 2)))
    timestep1 = 0.00001
    timestep2 = 0.000001

    # Create grains at z=0 and x and y and random values around the sheet
    grains = np.column_stack(((np.random.rand(number_of_grains, 2) - 0.5) * 2 * 0.5 ** 0.5 * sheet_radius,np.zeros(number_of_grains)))
    grains_type = sim.add_grains(grains, grain_diameter, 0.01)

    # I want to create a lammps geometry that will hold the grains. In an ideal world, what I would do is make
    # a cylinder that is very short that is the radius of the sheet. The problem is then that there may be grains
    # that are touching 3 sides of the cylinder, and LAMMPS flips if that happens for some reason. So instead I
    # will make three cylinders, each which comprises one of the walls that we want to contain the grains in
    wall1 = sim.add_walls("cylinder z 0 0 0.5 -0.002 0.2", type=grains_type, youngs_modulus=10 ** 3)
    wall2 = sim.add_walls("cylinder z 0 0 0.5 -0.2 0.002", type=grains_type, youngs_modulus=10 ** 3)
    wall3 = sim.add_walls(f"cylinder z 0 0 {sheet_radius} -0.2 0.2",type=grains_type,youngs_modulus=10 ** 3)

    # Make the grains repel each other
    sim.turn_on_granular_potential(grains_type, grains_type, 10 ** 3)

    # Since we aren't using the hardcore granular potential, there is nothing to slow down the grains,
    # and so we need some viscosity to chill them out
    grain_viscosity = sim.add_viscosity(0.001, grains_type)
    sim.design_dump_files(0.05, ["fx", "fy", "fz", "pressure", "bending_energy", "stretching_energy"])
    sim.run_simulation(0.5, timestep1)

    ## Now we will remove the walls and add the sheet
    sim.remove_something("walls",wall1)
    sim.remove_something("walls",wall2)
    sim.remove_something("walls",wall3)

    sim.remove_something('viscosity', grain_viscosity)

    sheet_type, sheet_bonds, sheet_angles = sim.add_circular_sheet(np.array([0,0,0]), sheet_radius,0.004, 0.002, 10**4, 0.5)
    sim.add_viscosity(0.00001, sheet_type)

    # I'm gonna reset the timestep manually here. The reason for this is that I really want to reset the timestep,
    # But also, if I wait until after I do the sim.move(), lammps is going to yell at me (detailed in the run_simulation) description
    sim.manually_edit_timestep(timestep2)

    # Find all of the particles which are around the edge of the sheet and make them not move
    clamped_particles = sim._particles.index[((sim._particles['x_position'] ** 2 + sim._particles['y_position'] ** 2) ** 0.5 > sheet_radius - 0.4 * grain_diameter) & (sim._particles['type']==2)].tolist()
    sim.move(particles = clamped_particles, xvel = 0, yvel = 0, zvel = 0)

    sim.connect_particles_to_elastic(grains_type, sheet_type, 1)

    sim.perturb(particles = [1,2], zdir = -1)
    # Here I don't have to worry about the timestep because I have run the sim once and called move(), so the timestep is ignored
    # Sorry, I know this is confusing, check out the description of run_simulation() to learn more
    sim.run_simulation(0.5)
