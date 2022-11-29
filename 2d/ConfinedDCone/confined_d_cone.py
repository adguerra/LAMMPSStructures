import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":

    # To use this file, enter the filename that you want the sim in in the command line python call, like "python confined_d_cone.py [filename]"
    simname = str(sys.argv[1])
    
    # Start a simulation with the name simname
    sim = Simulation(simname, 3, 0.18, 0.18, 0.08)

    # The 'Sheet thickness' will be what we use to set the energies (and therefore stretching and bending modulus)
    # but we'll use mesh particles of spacing 'mesh_particle_spacing' so the mesh isn't too small, and diameter
    # 'mesh_particle_diameter' so that the sheet is geometrically pretty thin but still slides against the cup.
    # You can change the spacing or diameter as you wish
    sheet_thickness = 1 * 0.127 * 10 ** -3
    mesh_particle_spacing = sheet_thickness * 10
    mesh_particle_diameter = sheet_thickness * 10

    Rc = 13 * 10 ** -3
    Rp = 55 * 10 ** -3
    Rr = 36 * 10 ** -3
    eps_max = 0.6
    E_sheet = 5 * 10 ** 5
    # We will offset the dcone from the center of the cup by this much, just for fun
    xoffset = 1 * sheet_thickness
    yoffset = 0 * sheet_thickness
    simtime = 1
    # These next three I got through a lot of trial and error. A larger timestep necessitates heavier, more viscous particles, but heavier 
    # and more viscous particles give the sheet momentum and viscosity, both of which work against the simulation being quasi-static
    density = 2 * 10 ** -2
    timestep = 5 * 10 ** -7
    viscosity = 1 * 10 ** -9

    # I want to create a lammps geometry that will be the cup. Problem is, I kind of want the lip of the
    # Cup to be a little rounded, which will allow the plate to slide. For that reason, I want to make
    # A for loop which adds a bunch of cylinders at different heights
    walls = []
    lip_thickness, lip_height = mesh_particle_spacing * 1.5, mesh_particle_spacing * 1
    cup_section_radii = np.linspace(Rr-lip_thickness / 2, Rr+lip_thickness / 2, 15).tolist()
    cup_section_heights = [-mesh_particle_diameter / 2 - lip_height * ((cpr-Rr)/(lip_thickness/2)) ** 2 for cpr in cup_section_radii]

    for rad,height in zip(cup_section_radii,cup_section_heights):
        walls.append(sim.add_walls(f"cylinder z 0 0 {rad} {-sim._height/2} {height} open 1 open 2", youngs_modulus=10 ** 4))

    # Create the sheet
    sheet_type, sheet_bonds, sheet_angles = sim.add_circular_sheet(
            np.array([xoffset,yoffset,0]),
            Rp,
            mesh_particle_spacing,
            sheet_thickness,
            E_sheet,
            density,
            mesh_particle_diameter = mesh_particle_diameter
    )

    # Clamp the clamp particles. Here we also clamp them in x and y, but this can be changed.
    clamp_particles = sim._particles.index[((sim._particles['x_position']-xoffset) ** 2 + (sim._particles['y_position']-yoffset) ** 2) ** 0.5 < Rc].tolist()
    sim.move(particles = clamp_particles, xvel = 0, yvel = 0, zvel = - eps_max * Rr / simtime)

    # Add the viscosity, which just helps the simulation from exploding, mimics the normal 
    # damping of air and slight viscoelasticity which we live in but don't appreciate
    sim.add_viscosity(viscosity, sheet_type)

    # Make the dump files and run the simulation
    sim.design_dump_files(0.005)
    sim.run_simulation(simtime, timestep)
