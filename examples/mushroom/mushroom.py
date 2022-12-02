import numpy as np
import pandas as pd
from lammpsWithPython.lammps_object import Simulation
import sys


if __name__ == "__main__":

    # To use this file, enter the filename that you want the sim in in the command line python call, like "python confined_d_cone.py [filename]"
    simname = str(sys.argv[1])

    n_beams = 10
    d_between_beams = 0.005
    beam_length = 0.095
    beam_thickness = 0.002
    n_particles_per_beam = 150
    squish_factor = 0.2
    E_beams = 0.96 * 10 ** 6
    E_walls = 10 ** 4

    simtime = 1
    density = 0.5
    viscosity = 2 * 10 ** -7
    timestep = 1 * 10 ** -7
    
    # Start a simulation with the name simname
    sim = Simulation(simname, 3, n_beams * d_between_beams, 0.01, 0.1)
    # Make the simulationStation hard. We can also do this sim periodically, so this is not required
    sim.add_walls(youngs_modulus = E_walls)

    # We're gonna do something tricky, which is we will turn on all granular interactions between all beams, and
    # Then in the loop where we add the beams we will turn off the granular interactions between the beams and themselves
    # If the particles in the beams are granular with themselves, it will mess up the beam mechanics
    sim.turn_on_granular_potential(youngs_modulus = E_walls)

    beam_positions = np.linspace(-0.5 * d_between_beams * (n_beams - 1), 0.5 * d_between_beams * (n_beams - 1), n_beams).tolist()
    for bp in beam_positions:
        beam, _, _ = sim.add_beam(n_particles_per_beam, np.array([bp,0,-beam_length/2]),np.array([bp,0,beam_length/2]),beam_thickness,E_beams,density)
        sim.turn_on_granular_potential(type1 = beam, type2 = beam, youngs_modulus = 0)

    # Clamp the clamp particles. Here we also clamp them in x and y, but this can be changed.
    bottom_clamp = [i for i in range(1,n_beams*n_particles_per_beam+1) if i%n_particles_per_beam in [1,2]]
    top_clamp = [i for i in range(1,n_beams*n_particles_per_beam+1) if i%n_particles_per_beam in [n_particles_per_beam-1,0]]
    sim.move(particles = bottom_clamp, xvel = 0, yvel = 0, zvel = squish_factor * beam_length / simtime)
    sim.move(particles = top_clamp, xvel = 0, yvel = 0, zvel = - squish_factor * beam_length / simtime)

    # Perturb the beams to buckle to the left or right randomly
    dirs = np.random.rand(len(beam_positions),1)
    p1 = sim.perturb(type = [i+1 for i in np.where(dirs>0.5)[0].tolist()],xdir = 1)
    p2 = sim.perturb(type = [i+1 for i in np.where(dirs<=0.5)[0].tolist()],xdir = -1)

    # Add the viscosity, which just helps the simulation from exploding, mimics the normal 
    # damping of air and slight viscoelasticity which we live in but don't appreciate
    sim.add_viscosity(viscosity)

    # Make the dump files and run the simulation
    sim.design_dump_files(0.01)
    sim.run_simulation(simtime, timestep)
