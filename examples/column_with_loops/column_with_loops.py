import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":
    simname = str(sys.argv[1])

    # Start a simulation with the name simname
    sim = Simulation(simname, 3, 0.4, 0.4, 0.4)

    grain_diameter = 0.013
    E_grains = 10 ** 6
    contact_properties = {
        "restitution" : 0.4,
        "poissons": 0.35,
        "xscaling": 1,
        "coeffric": 0.5,
        "gammar": 1,
        "rolfric": 0.5,
    }


    loop_thickness = 0.0015
    E_loop = 10 ** 4
    n_particles_per_loop = 180
    # This next thing can be used to independantly control the bending stiffness of the loops
    bending_multiplier = 1 

    loop_separation = 1.1 * grain_diameter
    number_of_grains = 375
    timestep = 0.0000005
    density = 0.01
    viscosity1 = 1 * 10 ** -6
    viscosity2 = 1 * 10 ** -9

    column_radius = 0.04
    column_height = 0.16

    # Create grains at x, y, and z at random points in a column with the radius of the final column and height the entire height of the simulation (so they settle)
    grains = np.column_stack((
            (np.random.rand(number_of_grains, 2) - 0.5) * 2 * 0.5 ** 0.5 * (column_radius - 0.9 * grain_diameter/2),
            (np.random.rand(number_of_grains) - 0.5) * (sim._height - grain_diameter/2)
            ))
    grains_type = sim.add_grains(grains, grain_diameter, density)

    # Make the simulation box and the column that the grains will settle in stiff
    sim.add_walls(youngs_modulus = E_grains, hardcore_dict = contact_properties)
    hold1 = sim.add_walls(f"cylinder z 0 0 {column_radius - loop_thickness} {-sim._height/2} {sim._height/2}", type=grains_type, youngs_modulus = E_grains)

    # Make the grains repel each other
    sim.turn_on_granular_potential(grains_type, grains_type, E_grains, contact_properties)
    sim.add_gravity()

    # Since we initialize some grains inside of one another, its nice to have viscosity at the beginning
    grain_viscosity = sim.add_viscosity(viscosity1, grains_type)
    sim.design_dump_files(0.05, ["pressure", "bending_energy", "stretching_energy"])
    sim.run_simulation(0.05, timestep)

    # Run it for another bit so the grains can settle
    sim.remove_something('viscosity', grain_viscosity)
    sim.run_simulation(0.45, timestep)

    ## Now we will remove the walls and add the loops, and add a wall outside of the loops
    sim.remove_something("walls", hold1)

    loop_types = []
    loop_heights = np.arange(-sim._height/2 + loop_separation, -sim._height/2 + column_height, loop_separation)
    for lh in loop_heights:
        lt,_,_ = sim.add_loop(n_particles_per_loop,np.array([0,0,lh]),column_radius-loop_thickness/2, np.array([0,0,1]),loop_thickness,E_loop,density,bending_multiplier=bending_multiplier)
        sim.turn_on_granular_potential(grains_type, lt,E_grains, contact_properties)
        loop_types.append(lt)
    
    
    hold2 = sim.add_walls(f"cylinder z 0 0 {column_radius} {-sim._height/2} {sim._height/2}", youngs_modulus = E_grains)
    sim.add_viscosity(viscosity2)
    sim.run_simulation(0.15, timestep)

    sim.remove_something("walls", hold2)
    sim.run_simulation(1, timestep)
