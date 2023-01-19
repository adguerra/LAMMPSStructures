import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":
    simname = str(sys.argv[1])

    sim = Simulation(simname, 3, 0.5, 0.5, 0.2)

    mesh_diameter = 0.004
    radius = mesh_diameter * 5
    height = 0.15
    timestep = 0.000001


    sheet_type, sheet_bonds, sheet_angles = sim.add_cylindrical_sheet(np.array([-0.075,0,0]), radius,height,np.array([1,0,0]),mesh_diameter, 0.002, 2 * 10 ** 4, 0.5)
    sim.add_viscosity(0.00001, type = sheet_type)

    # Create grains to squeeze the sheet
    n_grains = 20
    grain_d = 0.01
    grains1 = np.column_stack(((np.zeros(n_grains), np.linspace(-0.02,0.02,n_grains),np.zeros(n_grains) + grain_d / 2 + radius + mesh_diameter)))
    grains2 = np.column_stack(((np.zeros(n_grains), np.linspace(-0.02,0.02,n_grains),np.zeros(n_grains) - (grain_d / 2 + radius + mesh_diameter))))   

    grains_type1 = sim.add_grains(grains1, grain_d, 0.01)
    grains_type2 = sim.add_grains(grains2, grain_d, 0.01)

    sim.turn_on_granular_potential(sheet_type, youngs_modulus = 10 ** 3)
    move1 = sim.move(type = grains_type1,xvel = 0, yvel = 0,zvel = -radius)
    move2 = sim.move(type = grains_type2,xvel = 0, yvel = 0, zvel = radius)

    # sim.perturb(particles = [1,2], zdir = -1)
    sim.design_dump_files(0.025, ["fx", "fy", "fz", "pressure", "bending_energy", "stretching_energy"])
    sim.run_simulation(0.5, timestep)





