import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":
    simname = str(sys.argv[1])
    

    # Start a simulation with the name simname
    sim = Simulation(simname, 3, 0.5, 0.5, 0.2)

    grain_diameter = 0.025
    mesh_diameter = 0.004
    sheet_side_length = 0.159/1.5
    grain_packing_fraction = 0.85
    number_of_grains = int(np.round((sheet_side_length ** 2) * grain_packing_fraction / ((grain_diameter / 2) ** 2 * np.pi)))
    timestep1 = 0.00001
    timestep2 = 0.000001

    # Cohesive potential stuff
    grain_sheet_adhesion = 1

    # Create grains at z=0 and x and y and random values around the sheet
    grains = np.column_stack(((np.random.rand(number_of_grains, 2) - 0.5) * sheet_side_length,np.zeros(number_of_grains) + grain_diameter / 2 + 1.01 * mesh_diameter / 2))
    grains_type = sim.add_grains(grains, grain_diameter, 0.01)

    # I want to create a lammps geometry that will hold the grains. In an ideal world, what I would do is make
    # a box that is very short that is the size of the sheet. The problem is then that there may be grains
    # that are touching 3 sides of the box, and LAMMPS flips if that happens for some reason. So instead I
    # will make three boxes, each which comprises one of the walls in which we want to contain the grains
    wall1 = sim.add_walls("block -0.25 0.25 -0.25 0.25 -0.002 0.1", type=grains_type, youngs_modulus=10 ** 3)
    wall2 = sim.add_walls("block -0.25 0.25 -0.25 0.25 -0.1 0.002", type=grains_type, youngs_modulus=10 ** 3)
    wall3 = sim.add_walls(f"block -{sheet_side_length/2} {sheet_side_length/2} -{sheet_side_length/2} {sheet_side_length/2} -0.1 0.1",type=grains_type,youngs_modulus=10 ** 3)

    # Make the grains repel each other
    sim.turn_on_granular_potential(grains_type, grains_type, 10 ** 3)

    # Since we aren't using the hardcore granular potential, there is nothing to slow down the grains,
    # and so we need some viscosity to chill them out
    grain_viscosity = sim.add_viscosity(0.001, type = grains_type)
    sim.design_dump_files(0.05, ["fx", "fy", "fz", "pressure", "bending_energy", "stretching_energy"])
    sim.run_simulation(0.5, timestep1)

    ## Now we will remove the walls and add the sheet
    sim.remove_something("walls",wall1)
    sim.remove_something("walls",wall2)
    sim.remove_something("walls",wall3)

    sim.remove_something('viscosity', grain_viscosity)

    sheet_type, sheet_bonds, sheet_angles = sim.add_rectangular_sheet(np.array([0,0,0]), sheet_side_length,sheet_side_length,mesh_diameter, 0.002, 10 ** 4, 0.5)
    sim.add_viscosity(0.00001, type = sheet_type)

    
    sim.connect_particles_to_elastic(sheet_type, grain_sheet_adhesion, type=grains_type)

    # sim.perturb(particles = [1,2], zdir = -1)
    sim.run_simulation(0.5, timestep2)
