import os
from typing import List, Union, Tuple

import numpy as np
import pandas as pd


class Simulation:

    def __init__(
        self,
        simulation_name: str,
        dimension: int,
        length: float,
        width: float,
        height: float,
        x_bound: str = "f",
        y_bound: str = "f",
        z_bound: str = "f",
    ):
        """
        Here we are going to initialize a bunch of stuff about the simulation.
        Most of ths is random LAMMPS stuff that needs to be in order so that
        we can make an elastogranular sim.
        Inputs:
        - simulation_name: the name of the folder that this script will create and write files to
        - dimension: 2 or 3 for a 2d or 3d simulation
        - length, width, height: the length, width, and height of the box that we are going to run the simulation in
        - (x,y,z)_bound: whether the x, y, or z boundaries are fixed (f) periodic (p) or shrink-wrapped (s)
        """
        # TODO add 2d

        print(
            "\n### Making a simulation file in the folder '"
            + simulation_name
            + "' ###\n"
        )

        self._dimension = dimension
        self._length = length
        self._width = width
        self._height = height
        self._x_bound = x_bound
        self._y_bound = y_bound
        self._z_bound = z_bound
        self._simulation_name = simulation_name
        self._particles = pd.DataFrame(
            columns=[
                "type",
                "x_position",
                "y_position",
                "z_position",
                "diameter",
                "density",
            ]
        )
        self._num_sets_of_grains = 0
        self._num_beams = 0
        self._num_loops = 0
        self._num_sheets = 0
        self._type_iter = 0
        self._bond_type_iter = 0
        self._angle_type_iter = 0
        self._connections_iter = 0
        self._perturb_iter = 0
        self._visc_iter = 0
        self._group_iter = 0
        self._move_iter = 0
        self._wall_iter = 0
        self._walls = []
        self._timestep = False
        self._have_run = False
        self._dump_file_every_this_many_seconds = False

        # Make the directory that the lammps simulation will be placed in
        self._path = os.path.join(os.getcwd(), self._simulation_name)
        if not os.path.exists(self._path):
            os.mkdir(self._path)
        with open(os.path.join(self._path, "in.main_file"), "w") as f:
            # Atom style: angle, for bending and stretching rigidity, and sphere, for granular interactions
            f.write("atom_style hybrid angle sphere\n")
            # You can make this whatever you prefer! Consider lj for maximum eyebrow-raising potential
            f.write("units si\n")
            # We're just using one processor today
            f.write("processors	* * 1\n")
            # No idea what this does
            f.write("comm_modify vel yes\n")
            if dimension == 2:
                z_bound = "p"
                f.write("dimension 2")
            f.write("\n")
            # Define the box that the simulation is gonna happen in. First we create a region
            f.write(
                f"region simulationStation block -{length/2} {length/2} -{width/2} {width/2} -{height/2} {height/2}\n"
            )
            # Then we make it the simulation station. Here the numbers are:
            #   The number of types of atoms (100)
            #   The number of types of bonds (100000)
            #   The number of possible bonds per particle (100)
            #   The number of types of angles (10000)
            #   The number of possible angles per particle (100)
            # In this line we have to put something, and if we put too few, lammps will get mad at us. For example if we say we
            # are going to have 30 types of particles, but we insert 31 types if particles, LAMMPS will be like "thats more than we thought!"
            # and your sim will die. However, since we don't know how many we need at the time of the writing of this line,
            # We just put more that we think we'll need. And since we have lines like "bond_coeff * 0 0", we set all of the 
            # extra ones to zero. If you need to put in more than I allow here, feel free to change this line!
            f.write(
                "create_box 100 simulationStation bond/types 100000 extra/bond/per/atom 100 angle/types 10000 extra/angle/per/atom 100\n"
            )
            # Make the boundry fixed, periodic, or shrink wrapped
            f.write(
                "change_box	all boundary "
                + x_bound
                + " "
                + y_bound
                + " "
                + z_bound
                + "\n"
            )
            f.write("\n")
            # We will use a granular / lj pair style, giving us the ability to have granular interactions and also cohesive/repulsive energies
            f.write("pair_style hybrid/overlay granular lj/cut 0\n")
            # Set all of them to zero to start (LAMMPS doesnt like when anything is not set)
            f.write(
                "pair_coeff  * * granular hertz/material 0 0 0 tangential linear_nohistory 0 0\n"
            )
            f.write("pair_coeff  * * lj/cut 0 0 0\n")
            f.write("\n")
            # Not totally sure what this does
            f.write("special_bonds lj/coul 0 1.0 1.0\n")
            # These are the potentials we will use
            f.write("bond_style harmonic\n")
            f.write("angle_style cosine\n")
            f.write("\n")
            f.write("bond_coeff * 0 0\n")
            f.write("angle_coeff * 0\n")
            # Turn on integration
            f.write("fix integration all nve/sphere\n")
            # This allows you to put in more particles and bonds. You might need to change this if you want an even bigger simulation
            f.write("neigh_modify page 500000 one 50000\n")
            f.write("\n\n### Here we will begin to include particles ###\n\n")
        pass

    def add_grains(
        self, coords: np.array, diameter: float, density: float, filename: str = None
    ) -> int:
        """
        Add grains to the simulation
        Inputs:
        - coords: an Nx3 numpy array where N is the number of grains that we are trying to insert and
        the three columns are the x, y, and z coordinates of those grains
        - diameter: the diameter of the grains that you will be inserting
        - density: the density of the grains that you will be inserting
        - filename: the name of the file that you want to put the create lines in. I'd leave this empty
        Outputs:
        - The type id of the particles that you are adding. This is for input into methods like "add_pair_potential" or "remove_something"
        """

        self._type_iter += 1
        # Unless we are placing something specific, place grains
        if filename is None:
            print(f"\n### Adding {coords.shape[0]} grains to the simulation ###")
            self._num_sets_of_grains += 1
            filename = f"grains_{self._num_sets_of_grains}.txt"
        else:
            print(f"# Adding {coords.shape[0]} grains to the simulation #")
        # Add the "type" column to the coordinates as the first column of the array
        new_particles = np.vstack(
            (
                [self._type_iter] * coords.shape[0],
                coords.T,
                [diameter] * coords.shape[0],
                [density] * coords.shape[0],
            )
        ).T
        # Turn that into a dataframe
        if self._particles.shape[0]==0:
            indices = np.array(list(range(coords.shape[0])))+1
        else:
            indices = np.array(list(range(coords.shape[0])))+max(self._particles.index)+1
        new_particles = pd.DataFrame(new_particles, columns=self._particles.columns, index=indices)
        # Add those into the current set of particles
        self._particles = pd.concat([self._particles, new_particles], axis=0)
        # Put the particles into the simulation
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\ninclude " + filename + " \n")
            f.write(f"set type {self._type_iter} diameter {diameter}\n")
            f.write(f"set type {self._type_iter} density {density}\n")
        with open(os.path.join(self._path, filename), "w") as f:
            for i in range(coords.shape[0]):
                f.write(
                    f"create_atoms {self._type_iter} single {coords[i][0]} {coords[i][1]} {coords[i][2]}\n"
                )

        # Change the timestep to be able to handle these grains
        if self._timestep:
            self._timestep = min(self._timestep, 0.1 * (diameter / 9.8) ** 0.5)
        else:
            self._timestep = 0.01 * (diameter / 9.8) ** 0.5

        return self._type_iter

    def add_beam(
        self,
        n_particles: int,
        pos1: np.array,
        pos2: np.array,
        geometric_thickness: float,
        youngs_modulus: float,
        density: float,
        energetic_thickness_multiplier: float = 1,
        strecthing_multiplier: float = 1,
        bending_multiplier: float = 1,
    ) -> Tuple[int, int, int]:
        """
        Create a beam out of particles, replete with bonds and angles to give the beam mechanical properties
        Inputs:
        - n_particles: The number of particles that you want to make up your beam
        - pos1: A 1x3 np.array with the x, y, z coordinates of one end of the beam
        - pos2: A 1x3 np.array with the x, y, z coordinates of the other end of the beam
        - geometric_thickness: The diameter of the particles which make up the beam
        - youngs_modulus: The Young's modulus of the beam
        - density: The density of the particles which make up the beam. NOTE that this is different from
            the density of the beam itself
        - energetic_thickness_multiplier: This is to change the thickness term in the stretching and bending
            energies, which changes the stretching energy linearly and the bending energy cubically
        - stretching_multiplier: In case you want to manually edit the stretching energy
        - bending_multiplier: In case you want to manually edit the bending energy
        Outputs:
        - The number of the particles, bonds, and angles that make up this beam
        """
        if n_particles<3:
            raise Exception("Your beam must have more than 3 particles")

        print(f"\n### Creating a beam with {n_particles} particles ###")

        self._num_beams += 1
        # Make a list of the coordinates of the new particles
        coords = pos1
        for i in range(1, n_particles):
            coords = np.vstack((coords, pos1 + (pos2 - pos1) * i / (n_particles - 1)))
        # Make a list of the particles that we want to connect via bonds
        n_so_far = max(self._particles.index) if self._particles.index.size > 0 else 0
        tuplets = [[n_so_far + 1, n_so_far + 2]]
        for i in range(2, n_particles):
            tuplets.append([n_so_far + i, n_so_far + i + 1])
        tuplets = np.array(tuplets)
        # Make a list of the particles that we want to connect via angles
        triplets = [[n_so_far + 1, n_so_far + 2, n_so_far + 3]]
        for i in range(2, n_particles - 1):
            triplets.append([n_so_far + i, n_so_far + i + 1, n_so_far + i + 2])
        triplets = np.array(triplets)

        self.add_grains(
            coords, geometric_thickness, density, f"beam_{self._num_beams}.txt"
        )

        # Thickness
        h = geometric_thickness * energetic_thickness_multiplier
        # rest length of the bonds is the spacing between the particles
        rest_length = np.linalg.norm(pos2 - pos1) / (n_particles - 1)
        # dfactor is the diameter of a particle divided by the distance between particles in the beam
        dfactor = geometric_thickness / rest_length
        bond_stiffness = strecthing_multiplier * youngs_modulus * h * dfactor / 2
        angle_stiffness = bending_multiplier * youngs_modulus * (h ** 3) * dfactor / 12

        self.construct_many_bonds(tuplets, bond_stiffness, rest_length)
        self.construct_many_angles(triplets, angle_stiffness)

        if youngs_modulus > 0:
            # Change the timestep to be able to handle these grains
            m = (4 / 3) * np.pi * (geometric_thickness * 0.5) ** 3
            k = bond_stiffness * 2  # Somehow include angle_stiffness
            if self._timestep:
                self._timestep = min(self._timestep, ((1 / 2) * m / k) ** 0.5)  # TODO
            else:
                self._timestep = ((1 / 2) * m / k) ** 0.5

        return self._type_iter, self._bond_type_iter, self._angle_type_iter

    def add_loop(
        self,
        n_particles: int,
        center: np.array,
        radius: float,
        normal: np.array,
        geometric_thickness: float,
        youngs_modulus: float,
        density: float,
        energetic_thickness_multiplier: float = 1,
        strecthing_multiplier: float = 1,
        bending_multiplier: float = 1,
    ) -> Tuple[int, int, int]:
        """
        Create a loop out of particles, replete with bonds and angles to give the loop mechanical properties
        Inputs:
        - n_particles: The number of particles that you want to make up your beam
        - center: A 1x3 np.array with the x, y, z coordinates of the center of the loop
        - radius: The radius of the loop
        - normal: A 1x3 np.array which gives the vector normal to the plane which the loop occupies,
            for example, a ring which is laying flat on a table would have a normal \porm [0, 0, 1]
        - geometric_thickness: The diameter of the particles which make up the beam
        - youngs_modulus: The Young's modulus of the beam
        - density: The density of the particles which make up the beam. NOTE that this is different from
            the density of the beam itself
        - energetic_thickness_multiplier: This is to change the thickness term in the stretching and bending
            energies, which changes the stretching energy linearly and the bending energy cubically
        - stretching_multiplier: In case you want to manually edit the stretching energy
        - bending_multiplier: In case you want to manually edit the bending energy
        Outputs:
        - The number of the particles, bonds, and angles that make up this loop
        """
        self._num_beams += 1

        print(f"\n### Creating a loop with {n_particles} particles ###")

        # First I need to find two vectors that are perpendicular to the norm vector
        x1 = np.array([0, 0, 0])
        for i in range(3):
            if normal[i] != 0:
                x1[i - 1] = normal[i]
                x1[i] = -normal[i - 1]
                break
        x2 = np.cross(x1, normal)

        for vect in [normal, x1, x2]:
            vect = vect / np.linalg.norm(vect)

        # Now we start laying the coords
        coords = center + radius * x1
        thetas = np.linspace(2 * np.pi / n_particles, 2 * np.pi * (1 - 1 / n_particles), n_particles - 1)
        for theta in thetas:
            coords = np.vstack(
                (coords, center + radius * (x1 * np.cos(theta) + x2 * np.sin(theta)))
            )

        # Make a list of the particles that we want to connect via bonds
        n_so_far = max(self._particles.index) if self._particles.index.size > 0 else 0
        tuplets = [[n_so_far + n_particles, n_so_far + 1]]
        for i in range(1, n_particles):
            tuplets.append([n_so_far + i, n_so_far + i + 1])
        tuplets = np.array(tuplets)

        # Make a list of the particles that we want to connect via angles
        triplets = [[n_so_far + n_particles - 1, n_so_far + n_particles, n_so_far + 1]]
        triplets.append([n_so_far + n_particles, n_so_far + 1, n_so_far + 2])
        for i in range(1, n_particles - 1):
            triplets.append([n_so_far + i, n_so_far + i + 1, n_so_far + i + 2])
        triplets = np.array(triplets)

        self.add_grains(
            coords, geometric_thickness, density, f"beam_{self._num_beams}.txt"
        )

        # Thickness
        h = geometric_thickness * energetic_thickness_multiplier
        # rest length of the bonds is the spacing between the particles
        rest_length = np.linalg.norm(coords[0, :] - coords[1, :])
        # dfactor is the diameter of a particle divided by the distance between particles in the beam
        dfactor = geometric_thickness / rest_length
        bond_stiffness = strecthing_multiplier * youngs_modulus * h * dfactor / 2
        angle_stiffness = bending_multiplier * youngs_modulus * (h ** 3) * dfactor / 12
        self.construct_many_bonds(tuplets, bond_stiffness, rest_length)
        self.construct_many_angles(triplets, angle_stiffness)

        # Change the timestep to be able to handle these grains
        if youngs_modulus > 0:
            m = (4 / 3) * np.pi * (geometric_thickness * 0.5) ** 3
            k = bond_stiffness * 2  # Somehow include angle_stiffness
            if self._timestep:
                self._timestep = min(self._timestep, ((1 / 2) * m / k) ** 0.5)  # TODO
            else:
                self._timestep = ((1 / 2) * m / k) ** 0.5

        return self._type_iter, self._bond_type_iter, self._angle_type_iter

    def add_circular_sheet(
        self,
        center: np.array,
        radius: float,
        # normal: np.array,
        mesh_particle_spacing: float,
        energetic_thickness: float,
        youngs_modulus: float,
        density: float,
        mesh_particle_diameter: float = None,
        strecthing_multiplier: float = 1,
        bending_multiplier: float = 1,
    ) -> Tuple[int, int, int]:
        """
        Create a circular sheet out of particles, replete with bonds and angles to give the sheet mechanical properties
        NOTE: For now, all sheets have a normal of [0,0,1]. I might change this later (TODO)
        Inputs:
        - center: A 1x3 np.array with the x, y, z coordinates of the center of the sheet
        - radius: The radius of the sheet
        - mesh_particle_spacing: The resting distance between mesh particles in the sheet
        - energetic_thickness: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
            NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
            Explicitly, for you elasticity mathematicians out there, this is h.
        - youngs_modulus: The Young's modulus of the sheet material
        - density: The density of the particles which make up the sheet. NOTE that this is different from
            the density of the sheet itself
        - mesh_particle_diameter: This term changes the diameter of the mesh particles. If this is not set (left set to None),
            the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
            will look like they're right next to each other
        - stretching_multiplier: In case you want to manually edit the stretching energy
        - bending_multiplier: In case you want to manually edit the bending energy
        Outputs:
        - The number of the particles, bonds, and angles that make up this sheet
        """
        self._num_sheets += 1

        geometric_thickness = mesh_particle_diameter if mesh_particle_diameter else mesh_particle_spacing

        coords = [[0, 0, 0]]
        for i in range(
            round(-radius * 2 / mesh_particle_spacing), round(radius * 2 / mesh_particle_spacing)
        ):
            for j in range(
                round(-radius * 2 / mesh_particle_spacing),
                round(radius * 2 / mesh_particle_spacing),
            ):
                if i % 2 == 0:
                    x = i * np.cos(np.pi / 6) * mesh_particle_spacing
                    y = j * mesh_particle_spacing
                else:
                    x = i * np.cos(np.pi / 6) * mesh_particle_spacing
                    y = (j - np.sin(np.pi / 6)) * mesh_particle_spacing
                if (x ** 2 + y ** 2) ** 0.5 < radius:
                    coords.append(list(center + np.array([x, y, 0])))
        del coords[0]
        coords = np.array(coords)

        print(f"\n### Creating a circular sheet with {coords.shape[0]} particles ###")

        n_so_far = max(self._particles.index) if self._particles.index.size > 0 else 0
        tuplets = [[0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j)< mesh_particle_spacing * 1.01:
                    tuplets.append([n_so_far + i + 1, n_so_far + j + 1])
                    a += 1
                    if a == 3:break
        del tuplets[0]
        tuplets = np.array(tuplets)

        triplets = [[0, 0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j) < mesh_particle_spacing * 1.01:
                    for k in range(j + 1, coords.shape[0]):
                        if (
                            _dist(coords,j,k) < mesh_particle_spacing * 1.01
                            and _dist(coords,i,k) > mesh_particle_spacing * 1.9
                        ):
                            triplets.append([n_so_far + i + 1, n_so_far + j + 1, n_so_far + k + 1])
                            break
                    a += 1
                    if a == 3:break

        del triplets[0]
        triplets = np.array(triplets)

        self.add_grains(
            coords, geometric_thickness, density, f"sheet_{self._num_sheets}.txt"
        )
        # Thickness
        h = energetic_thickness
        # rest length of the bonds is the spacing between the particles
        rest_length = mesh_particle_spacing
        # dfactor is the diameter of a particle divided by the distance between particles in the beam
        bond_stiffness = (
            strecthing_multiplier * youngs_modulus * h * (3 ** 0.5) / 4
        )
        angle_stiffness = (
            (4 / (3 * (3 ** 2)))
            * bending_multiplier
            * youngs_modulus
            * (h ** 3)
            / (12 * (1 - (1 / 3) ** 2))
        )

        self.construct_many_bonds(tuplets, bond_stiffness, rest_length)
        self.construct_many_angles(triplets, angle_stiffness)

        if youngs_modulus > 0:
            m = (4 / 3) * np.pi * (geometric_thickness * 0.5) ** 3
            k = bond_stiffness * 2  # Somehow include angle_stiffness
            if self._timestep:
                self._timestep = min(self._timestep, ((1 / 6) * m / k) ** 0.5)  # TODO
            else:
                self._timestep = ((1 / 6) * m / k) ** 0.5

        return self._type_iter, self._bond_type_iter, self._angle_type_iter

    def add_rectangular_sheet(
        self,
        center: np.array,
        side_length1: float,
        side_length2: float,
        # normal: np.array,
        mesh_particle_spacing: float,
        energetic_thickness: float,
        youngs_modulus: float,
        density: float,
        mesh_particle_diameter: float = None,
        strecthing_multiplier: float = 1,
        bending_multiplier: float = 1,
    ) -> Tuple[int, int, int]:
        """
        Create a rectangular sheet out of particles, replete with bonds and angles to give the sheet mechanical properties
        NOTE: For now, all sheets have a normal of [0,0,1]. I might change this later (TODO)
        Inputs:
        - center: A 1x3 np.array with the x, y, z coordinates of the center of the sheet
        - side_length1: The x-dimension length of the sheet
        - side_length2: The y-dimension length of the sheet
        - mesh_particle_spacing: The resting distance between mesh particles in the sheet
        - energetic_thickness: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
            NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
            Explicitly, for you elasticity mathematicians out there, this is h.
        - youngs_modulus: The Young's modulus of the sheet material
        - density: The density of the particles which make up the sheet. NOTE that this is different from
            the density of the sheet itself
        - mesh_particle_diameter: This term changes the diameter of the mesh particles. If this is not set (left set to None),
            the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
            will look like they're right next to each other
        - stretching_multiplier: In case you want to manually edit the stretching energy
        - bending_multiplier: In case you want to manually edit the bending energy
        Outputs:
        - The number of the particles, bonds, and angles that make up this sheet
        """
        self._num_sheets += 1

        geometric_thickness = mesh_particle_diameter if mesh_particle_diameter else mesh_particle_spacing

        coords = [[0, 0, 0]]
        for i in range(
            round(- side_length1 * 1.5 / mesh_particle_spacing), round(side_length1 * 1.5 / mesh_particle_spacing)
        ):
            for j in range(
                round(-side_length2 * 1.5 / mesh_particle_spacing),
                round(side_length2 * 1.5 / mesh_particle_spacing),
            ):
                if i % 2 == 0:
                    x = i * np.cos(np.pi / 6) * mesh_particle_spacing
                    y = j * mesh_particle_spacing
                else:
                    x = i * np.cos(np.pi / 6) * mesh_particle_spacing
                    y = (j - np.sin(np.pi / 6)) * mesh_particle_spacing
                if abs(x) < side_length1 / 2 and abs(y) < side_length2 / 2 :
                    coords.append(list(center + np.array([x, y, 0])))
        del coords[0]
        coords = np.array(coords)

        print(f"\n### Creating a circular sheet with {coords.shape[0]} particles ###")

        n_so_far = max(self._particles.index) if self._particles.index.size > 0 else 0
        tuplets = [[0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j)< mesh_particle_spacing * 1.01:
                    tuplets.append([n_so_far + i + 1, n_so_far + j + 1])
                    a += 1
                    if a == 3:break
        del tuplets[0]
        tuplets = np.array(tuplets)

        triplets = [[0, 0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j) < mesh_particle_spacing * 1.01:
                    for k in range(j + 1, coords.shape[0]):
                        if (
                            _dist(coords,j,k) < mesh_particle_spacing * 1.01
                            and _dist(coords,i,k) > mesh_particle_spacing * 1.9
                        ):
                            triplets.append([n_so_far + i + 1, n_so_far + j + 1, n_so_far + k + 1])
                            break
                    a += 1
                    if a == 3:break

        del triplets[0]
        triplets = np.array(triplets)

        self.add_grains(
            coords, geometric_thickness, density, f"sheet_{self._num_sheets}.txt"
        )
        # Thickness
        h = energetic_thickness
        # rest length of the bonds is the spacing between the particles
        rest_length = mesh_particle_spacing
        # dfactor is the diameter of a particle divided by the distance between particles in the beam
        bond_stiffness = (
            strecthing_multiplier * youngs_modulus * h * (3 ** 0.5) / 4
        )
        angle_stiffness = (
            (4 / (3 * (3 ** 2)))
            * bending_multiplier
            * youngs_modulus
            * (h ** 3)
            / (12 * (1 - (1 / 3) ** 2))
        )

        self.construct_many_bonds(tuplets, bond_stiffness, rest_length)
        self.construct_many_angles(triplets, angle_stiffness)

        if youngs_modulus > 0:
            m = (4 / 3) * np.pi * (geometric_thickness * 0.5) ** 3
            k = bond_stiffness * 2  # Somehow include angle_stiffness
            if self._timestep:
                self._timestep = min(self._timestep, ((1 / 6) * m / k) ** 0.5)  # TODO
            else:
                self._timestep = ((1 / 6) * m / k) ** 0.5

        return self._type_iter, self._bond_type_iter, self._angle_type_iter

    def add_cylindrical_sheet(
        self,
        bottom_center: np.array,
        radius: float,
        height: float,
        normal: np.array,
        mesh_particle_spacing: float,
        energetic_thickness: float,
        youngs_modulus: float,
        density: float,
        mesh_particle_diameter: float = None,
        strecthing_multiplier: float = 1,
        bending_multiplier: float = 1,
    ) -> Tuple[int, int, int]:
        """
        Create a cylindrical sheet out of particles, replete with bonds and angles to give the sheet mechanical properties
        Inputs:
        - center: A 1x3 np.array with the x, y, z coordinates of the center of the bottom ring of the cylinder
        - radius: The radius of the sheet
        - normal: The vector which points down the center of the cylinder
        - mesh_particle_spacing: The resting distance between mesh particles in the sheet
        - energetic_thickness: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
            NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
            Explicitly, for you elasticity mathematicians out there, this is h.
        - youngs_modulus: The Young's modulus of the sheet material
        - density: The density of the particles which make up the sheet. NOTE that this is different from
            the density of the sheet itself
        - mesh_particle_diameter: This term changes the diameter of the mesh particles. If this is not set (left set to None),
            the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
            will look like they're right next to each other
        - stretching_multiplier: In case you want to manually edit the stretching energy
        - bending_multiplier: In case you want to manually edit the bending energy
        Outputs:
        - The number of the particles, bonds, and angles that make up this sheet
        """
        self._num_sheets += 1

        geometric_thickness = mesh_particle_diameter if mesh_particle_diameter else mesh_particle_spacing

        # First I need to find two vectors that are perpendicular to the norm vector
        x1 = np.array([0, 0, 0])
        for i in range(3):
            if normal[i] != 0:
                x1[i - 1] = normal[i]
                x1[i] = -normal[i - 1]
                break
        x2 = np.cross(x1, normal)

        for vect in [normal, x1, x2]:
            vect = vect / np.linalg.norm(vect)

        dth = 2 * np.pi / np.round(2 * np.pi * radius / (mesh_particle_spacing * 3 ** 0.5))
        thetas1 = np.arange(0,2 * np.pi,dth)
        thetas2 = thetas1 + dth / 2
        x = radius * np.cos(thetas1)
        y = radius * np.sin(thetas1)
        actuald = (((x[1]-x[0]) ** 2 + (y[1]-y[0]) **2) ** 0.5)/(3 ** 0.5)
        print(actuald)

        numrings = round(2*height/actuald)
        coords = [[0,0,0]]
        for ring in range(numrings):
            thetas = thetas1 if ring % 2 == 0 else thetas2
            for theta in thetas:
                coords.append(list(bottom_center + radius * (x1 * np.cos(theta) + x2 * np.sin(theta)) + normal * ring * actuald / 2))
        del coords[0]
        coords = np.array(coords)

        print(f"\n### Creating a circular sheet with {coords.shape[0]} particles ###")

        mesh_particle_spacing = actuald

        n_so_far = max(self._particles.index) if self._particles.index.size > 0 else 0
        tuplets = [[0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j)< mesh_particle_spacing * 1.01:
                    tuplets.append([n_so_far + i + 1, n_so_far + j + 1])
                    a += 1
                    if a == 3:break
        del tuplets[0]
        tuplets = np.array(tuplets)

        triplets = [[0, 0, 0]]
        for i in range(coords.shape[0]):
            a = 0
            for j in range(i + 1, coords.shape[0]):
                if _dist(coords,i,j) < mesh_particle_spacing * 1.05:
                    for k in range(j + 1, coords.shape[0]):
                        if (
                            _dist(coords,j,k) < mesh_particle_spacing * 1.05
                            and _dist(coords,i,k) > mesh_particle_spacing * 1.9
                        ):
                            triplets.append([n_so_far + i + 1, n_so_far + j + 1, n_so_far + k + 1])
                            break
                    a += 1
                    if a == 3:break

        del triplets[0]
        triplets = np.array(triplets)

        self.add_grains(
            coords, geometric_thickness, density, f"sheet_{self._num_sheets}.txt"
        )
        # Thickness
        h = energetic_thickness
        # rest length of the bonds is the spacing between the particles
        rest_length = mesh_particle_spacing
        # dfactor is the diameter of a particle divided by the distance between particles in the beam
        bond_stiffness = (
            strecthing_multiplier * youngs_modulus * h * (3 ** 0.5) / 4
        )
        angle_stiffness = (
            (4 / (3 * (3 ** 2)))
            * bending_multiplier
            * youngs_modulus
            * (h ** 3)
            / (12 * (1 - (1 / 3) ** 2))
        )

        self.construct_many_bonds(tuplets, bond_stiffness, rest_length)
        self.construct_many_angles(triplets, angle_stiffness)

        if youngs_modulus > 0:
            m = (4 / 3) * np.pi * (geometric_thickness * 0.5) ** 3
            k = bond_stiffness * 2  # Somehow include angle_stiffness
            if self._timestep:
                self._timestep = min(self._timestep, ((1 / 6) * m / k) ** 0.5)  # TODO
            else:
                self._timestep = ((1 / 6) * m / k) ** 0.5

        return self._type_iter, self._bond_type_iter, self._angle_type_iter

    def construct_many_bonds(
        self, tuplets: np.array, stiffness: float, rest_length: float
    ) -> int:
        """
        Add harmonic bonds between particles.
        Inputs:
        - tuplets: An Nx2 np.array where the rows are pairs of particles that you want to make a bond between
        - stiffness: The stiffness of these bonds
        - rest_length: The rest length of the bonds
        Outputs:
        - The number of these bonds
        """
        self._bond_type_iter += 1

        print(f"# Creating {tuplets.shape[0]} Bonds #")

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"bond_coeff {self._bond_type_iter} {stiffness} {rest_length}\n")
            f.write(f"include bonds_{self._bond_type_iter}.txt\n")
        with open(
            os.path.join(self._path, f"bonds_{self._bond_type_iter}.txt"), "w"
        ) as f:
            for i in range(tuplets.shape[0]):
                f.write(
                    f"create_bonds single/bond {self._bond_type_iter} {tuplets[i][0]} {tuplets[i][1]}\n"
                )

        return self._bond_type_iter

    def construct_many_angles(self, triplets: np.array, stiffness: float) -> int:
        """
        Add cosine angles between particles.
        Inputs:
        - triplets: An Nx3 np.array where the rows are triplets of particles that you want to make an angle between
        - stiffness: The stiffness of these angles
        Outputs:
        - The number of these angles
        """
        self._angle_type_iter += 1

        print(f"# Creating {triplets.shape[0]} Angles #")

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"angle_coeff {self._angle_type_iter} {stiffness}\n")
            f.write(f"include angles_{self._angle_type_iter}.txt\n")
        with open(
            os.path.join(self._path, f"angles_{self._angle_type_iter}.txt"), "w"
        ) as f:
            for i in range(triplets.shape[0]):
                f.write(
                    f"create_bonds single/angle {self._angle_type_iter} {triplets[i][0]} {triplets[i][1]} {triplets[i][2]}\n"
                )

        return self._angle_type_iter

    def turn_on_granular_potential(
        self, type1: int = None, type2: int = None, youngs_modulus: float = None, hardcore_dict: dict = None
    ):
        # restitution = None, poissons = None, xscaling = None, coeffric = None
        """
        Make two types of particles interact in a granular way. This can either be a simple
        contact potential, which repels particles which are overlapping, or it can be a super
        complicated potential which adds all sorts of granular mechanics. If you want the chill
        potentiall, all you need to input is:
        - type1, type2: The types of particles you want to add an interaction potential to. This is output from methods like "add_grains," "add_beam," etc. 
        - youngs_modulus: The youngs modulus of the particles -- this will determine how much they can overlap given some confining force

        If you don't pass in type2, it will turn on the granular potential between type1 and all other particles. If you
        pass in neither type1 nor type2 it will turn on the granular potential between all particles.
        
        If you want to add additional physics, such that these particles actually
        behave like some sort of granular material, then you should:
        - input a "hardcore_dict" which contains:
            - restitution: The restitution coefficient
            - poissons: The Poisson's ratio of the interaction
            - xscaling: A scalar multiplier for the tangential damping
            - coeffric: The coefficient of sliding friction
            - gammar: The degree of rolling damping
            - rolfric: The coefficient of rolling friction
        I'd check out the lammps documentation for the "pair_style granular command" if you wanna be hardcore. Note that the downside
        to being hardcore is that it makes the simulation take much longer
        """

        if hardcore_dict and not all(
            key in hardcore_dict
            for key in (
                "restitution",
                "poissons",
                "xscaling",
                "coeffric",
                "gammar",
                "rolfric",
            )
        ):
            raise Exception(
                "If the hardcore flag is on, you also need to input 'restitution', 'poissons', 'xscaling', 'coeffric', 'gammar', 'rolfric' in the hardcore_dict"
            )

        if not type1 is None and not type2 is None:
            a = [type1, type2]
            type1 = min(a)
            type2 = max(a)

        print(
            f"\n### Initiating a granular potential between type {'all' if type1 is None else type1} and type {'all' if type2 is None else type2} Particles ###"
        )

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            if hardcore_dict:
                kr = (
                    4
                    * youngs_modulus
                    / (
                        2
                        * (2 - hardcore_dict["poissons"])
                        * (1.0 + hardcore_dict["poissons"])
                    )
                )
                f.write(f"pair_coeff {'*' if type1 is None else type1} {'*' if type2 is None else type2} granular &\n")
                f.write(
                    f"hertz/material {youngs_modulus} {hardcore_dict['restitution']} {hardcore_dict['poissons']} tangential mindlin NULL {hardcore_dict['xscaling']} {hardcore_dict['coeffric']} &\n"
                )
                f.write(
                    f"rolling sds {kr} {hardcore_dict['gammar']} {hardcore_dict['rolfric']} twisting marshall damping tsuji\n"
                )
            else:
                f.write(
                    f"pair_coeff {'*' if type1 is None else type1} {'*' if type2 is None else type2} granular hertz/material {youngs_modulus} 0 0.5 tangential linear_nohistory 0 0\n"
                )
        pass

    def turn_on_cohesive_potential(
        self, type1: int = None, type2: int = None, cohesive_strength: float = None, rest_length: float = None, cutoff:float = None
    ):
        """
        Make two types of particles have a cohesive (or repellant) lj potential. This takes in:
        - type1, type2: The types of particles you want to add an interaction potential to. This is output from methods like "add_grains," "add_beam," etc. 
        - cohesive_strength: This is epsilon in https://docs.lammps.org/pair_lj.html
        - rest_length: The rest length of the bonds. This is NOT sigma in https://docs.lammps.org/pair_lj.html
        - cutoff: If two particles get farther than "cutoff" apart, they will not be coherent any more. This can just be really big,
            or like 2.5 times rest_length, where the potential would pretty much disappear anyways

        If you don't pass in type2, it will turn on the potential between type1 and all other particles. If you
        pass in neither type1 nor type2 it will turn on the potential between all particles.
        
        """

        if any(thing is None for thing in [cohesive_strength, rest_length, cutoff]):
            raise Exception(
                "You must pass in a bond_strength, rest_length, and cutoff, sorry :/. The cutoff can just be really big if you don't actually want one"
            )

        if not type1 is None and not type2 is None:
            a = [type1, type2]
            type1 = min(a)
            type2 = max(a)

        print(
            f"\n### Initiating an lj potential between type {'all' if type1 is None else type1} and type {'all' if type2 is None else type2} Particles ###"
        )

        sigma = rest_length / ( 2 ** (1/6) )

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(
                f"pair_coeff {'*' if type1 is None else type1} {'*' if type2 is None else type2} lj/cut {cohesive_strength} {sigma} {cutoff}\n"
            )
        pass

    def connect_particles_to_elastic(
        self,
        elastic_type: int,
        stiffness: float,
        particles: List[int] = None,
        type: int = None,
        n_bonds_per_grain: int = 3,
        cutoff: float = None,
    ) -> int:
        """
        TODO: Add cutoff

        Create bonds between a set of particles and an elastic material. This will set the rest length of the bonds equal to the
        distance, at the time of the command, between each particle and the surface.
        Inputs:
        EITHER
        - particles: a list of particles
        OR 
        - type: a type of particles
        and:
        - elastic_type: The type of the elastic
        - stiffness: The stiffness of the bonds
        - n_bonds_per_particle: How many connections you want between each particle and the elastic. This is to restrict degrees of freedom. This can be 1 2 or 3.

        Note that this cannot be undone once it is done, unless you delete the particles or the surface
        """
        self._connections_iter += 1

        if n_bonds_per_grain > 3:
            raise Exception(f"Can't do more than 3 bonds per grain, you tried to do {n_bonds_per_grain}")
        
        if not particles is None and not type is None:
            raise Exception("Can either pass in a list of particles, or a type of particles, not both")

        ## If we get a type of particles, turn particles into a list of the particles of that type
        if not type is None:
            print(
                f"\n### Connecting particles of type {type} to the surface of type {elastic_type} ###"
            )
            particles = self._particles.index[
                self._particles["type"] == type
            ].tolist()
        else:
            print(
                f"\n### Connecting A list of particles to the surface of type {elastic_type} ###"
            )
        #print(particles)
        filename = f"connect_{self._connections_iter}.txt"
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\ninclude " + filename + " \n")

        with open(os.path.join(self._path, filename), "w") as f:
            f.write(f"group elastic_connect type {elastic_type}")
            for particle in particles:
                f.write(f"\n\n\ngroup now id {particle}")
                f.write("\nvariable distances atom abs(sqrt((xcm(now,x)-x)^2+(xcm(now,y)-y)^2+(xcm(now,z)-z)^2))")
                f.write("\nvariable rnow equal sqrt((xcm(now,x))^2+(xcm(now,y))^2)")
                f.write("\ncompute cmindist elastic_connect reduce min v_distances")
                f.write("\nthermo_style	custom step atoms time c_cmindist")
                f.write("\nrun 0")
                f.write("\nthermo_style custom step atoms time")
                f.write("\nvariable mindist equal c_cmindist")
                f.write("\nvariable rhinow equal ${mindist}*1.0000001")
                self._bond_type_iter += 1
                f.write(f"\nbond_coeff {self._bond_type_iter} {stiffness} " + "${mindist}")
                f.write(f"\ncreate_bonds many now elastic_connect {self._bond_type_iter} 0 "+"${rhinow}")
                if n_bonds_per_grain > 1:
                    f.write("\nvariable xone equal xcm(now,x)")
                    f.write("\nvariable yone equal xcm(now,y)")
                    f.write("\nvariable zone equal xcm(now,z)")
                    f.write("\nregion roneplus sphere ${xone} ${yone} ${zone} ${rhinow}")
                    f.write("\ngroup goneplus region roneplus")
                    f.write("\ngroup notTheClosest subtract all goneplus")
                    f.write("\nuncompute cmindist")
                    f.write("\ncompute cmindist notTheClosest reduce min v_distances")
                    f.write("\nthermo_style custom step atoms time c_cmindist")
                    f.write("\nrun 0")
                    f.write("\nthermo_style custom step atoms time")
                    f.write("\nvariable mindist equal c_cmindist")
                    f.write("\nvariable rhinow equal ${mindist}*1.0000001")
                    self._bond_type_iter += 1
                    f.write(f"\nbond_coeff {self._bond_type_iter} {stiffness} " + "${mindist}")
                    f.write(f"\ncreate_bonds many now notTheClosest {self._bond_type_iter} 0 "+"${rhinow}")
                    if n_bonds_per_grain > 2:
                        f.write("\nregion rtwoplus sphere ${xone} ${yone} ${zone} ${rhinow}")
                        f.write("\ngroup gtwoplus region rtwoplus")
                        f.write("\ngroup notTheClosestOrTheSecondClosest subtract all gtwoplus")
                        f.write("\nuncompute cmindist")
                        f.write("\ncompute cmindist notTheClosestOrTheSecondClosest reduce min v_distances")
                        f.write("\nthermo_style custom step atoms time c_cmindist")
                        f.write("\nrun 0")
                        f.write("\nthermo_style custom step atoms time")
                        f.write("\nvariable mindist equal c_cmindist")
                        f.write("\nvariable rhinow equal ${mindist}*1.0000001")
                        self._bond_type_iter += 1
                        f.write(f"\nbond_coeff {self._bond_type_iter} {stiffness} " + "${mindist}")
                        f.write(f"\ncreate_bonds many now notTheClosestOrTheSecondClosest {self._bond_type_iter} 0 "+"${rhinow}")
                f.write("\nuncompute cmindist")
                f.write("\ngroup now delete")
                if n_bonds_per_grain > 2:
                    f.write("\ngroup notTheClosestOrTheSecondClosest delete")
                    f.write("\ngroup gtwoplus delete")
                    f.write("\nregion rtwoplus delete")
                if n_bonds_per_grain > 1:
                    f.write("\ngroup notTheClosest delete")
                    f.write("\ngroup goneplus delete")
                    f.write("\nregion roneplus delete")
            f.write(f"\ngroup elastic_connect delete")
        return self._connections_iter

    def add_walls(
        self,
        region_details: str = None,
        particles: Union[List[int],int] = None,
        type: Union[List[int],int] = None,
        youngs_modulus: float = None,
        hardcore_dict: dict = None,
    ) -> int:
        """
        Make a lammps "region", and optionally make the walls of that region have a granular potential with a type of particles.
        This can either be a simple contact potential, which repels particles which are overlapping, or it can be a super
        complicated potential which adds all sorts of granular mechanics. To make the region, you have to add the
        - region_details: Everything that is required to define a lammps region other than the region name. 
        For example, if I wanted to define a cylinder which was aligned with the z axis and whose center was at x=0, y=0
        with a radius of 0.5 and which went from z=-0.2 to z=0.2 then I would have region_details = "cylinder z 0 0 0.5 -0.2 0.2".
        The full set of options are detailed here: https://docs.lammps.org/region.html
        
        If you pass in a youngs_modulus, the walls of the region that you are defining will become stiff to 
        some or all of the grains. Further options are:
        - If you pass in neither a particle or list of particle ids (particles) nor a type or list of types of particles,
            ALL particles in the simulation will interact with the walls
        - If you pass in either particles, or type, the particles or type of particles that you pass in will interact with the walls
        - If you want to add additional physics, such that these particles actually behave like some sort of granular material, then you should input a "hardcore_dict" which contains:
            - restitution: The restitution coefficient
            - poissons: The Poisson's ratio of the interaction
            - xscaling: A scalar multiplier for the tangential damping
            - coeffric: The coefficient of sliding friction
            - gammar: The degree of rolling damping
            - rolfric: The coefficient of rolling friction
        I'd check out the lammps documentation for the "pair_style granular command" if you wanna be hardcore. Note that the downside
        to being hardcore is that it makes the simulation take much longer

        One more thing -- if you pass in nothing for the region details, this will make the walls of the simulation station hard to the
        particles or type that you input, with youngs modulus and hardcore dict acting as usual. Note that this cannot currently be undone (TODO)
        """
        if hardcore_dict and not all(
            key in hardcore_dict
            for key in (
                "restitution",
                "poissons",
                "xscaling",
                "coeffric",
                "gammar",
                "rolfric",
            )
        ):
            raise Exception(
                "If the hardcore flag is on, you also need to input 'restitution', 'poissons', 'xscaling', 'coeffric', 'gammar', 'rolfric' in the hardcore_dict"
            )

        if not particles is None and not type is None:
            raise Exception(
                "You can make EITHER a particle OR a type of particle OR all of the particles in the simulation (input neither particles nor type) interact with a wall"
            )
        
        if not type is None:
            if isinstance(type,int):
                group = f"type {type}"
            else: 
                group = f"type " + " ".join([str(typ) for typ in type])
        elif not particles is None:
            if isinstance(particles,int):
                group = f"id {particles}"
            else: 
                group = f"id " + " ".join([str(part) for part in particles])
        else: group = "all"

        self._wall_iter += 1
        self._walls.append([bool(particles or type), bool(youngs_modulus)])

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write("\n")
            if group != "all":
                f.write(f"group group_wall_{self._wall_iter} " + group + "\n")
                group = f"group_wall_{self._wall_iter}"

            if not region_details is None:
                f.write(f"region region_{self._wall_iter} " + region_details + "\n")
                region = f"region_{self._wall_iter}"
            else:
                region = "simulationStation"

            if youngs_modulus:
                if hardcore_dict:
                    kr = (
                        4
                        * youngs_modulus
                        / (
                            2
                            * (2 - hardcore_dict["poissons"])
                            * (1.0 + hardcore_dict["poissons"])
                        )
                    )
                    f.write(
                        f"fix fix_walls_{self._wall_iter} "
                        + group
                        + " wall/gran/region granular &\n"
                    )
                    f.write(
                        f"hertz/material {youngs_modulus} {hardcore_dict['restitution']} {hardcore_dict['poissons']} tangential mindlin NULL {hardcore_dict['xscaling']} {hardcore_dict['coeffric']} &\n"
                    )
                    f.write(
                        f"rolling sds {kr} {hardcore_dict['gammar']} {hardcore_dict['rolfric']} twisting marshall damping tsuji region " + region + "\n"
                    )
                else:
                    f.write(
                        f"fix fix_walls_{self._wall_iter} "
                        + group
                        + f" wall/gran/region granular hertz/material {youngs_modulus} 0.25 0.25 tangential linear_nohistory 0 0 region " + region + "\n"
                    )
        return self._wall_iter

    def move(
        self,
        particles: Union[List[int], int] = None,
        type: Union[List[int],int] = None,
        xvel: float = None,
        yvel: float = None,
        zvel: float = None,
    ) -> int:
        """
        Move a set of particles at some velocity
        Pass in EITHER:
        - particles: A partice or list of particles to move
        OR
        - type: A type or list of types of particle to move (output from methods like "add_grains")

        And:
        xvel, yvel, zvel: The velocity in the x, y, and z directions to move that particle

        Note: If you pass in 'None' to either xvel, yvel, or zvel, or leave them blank, those 
        velocities will not be mandated, and the particles will be free to move in those directions
        """
        if not particles is None and not type is None:
            raise Exception("You can move EITHER a particle OR a type of particle")
        if particles is None and type is None:
            raise Exception("You must move either a particle or a type of particle")
        
        if not type is None:
            if isinstance(type,int):
                group = f"type {type}"
            else: 
                group = f"type " + " ".join([str(typ) for typ in type])
        else:
            if isinstance(particles,int):
                group = f"id {particles}"
            else: 
                group = f"id " + " ".join([str(part) for part in particles])

        self._move_iter += 1

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\ngroup group_move_{self._move_iter} " + group + "\n")
            f.write(
                f"fix fix_move_{self._move_iter} group_move_{self._move_iter} move linear "
                + " ".join(str(dir_vel) if not dir_vel is None else 'NULL' for dir_vel in [xvel, yvel, zvel])
                + "\n"
            )

        return self._move_iter

    def add_gravity(self, magnitude=9.8, xdir=0, ydir=0, zdir=-1):
        """
        Add gravity to the simulation
        Inputs:
        - magnitude: magnitude of the gravitational force
        - xdir, ydir, zdir: direction of the gravitational force
        """
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\nfix grav all gravity {magnitude} vector {xdir} {ydir} {zdir}\n")
        pass

    def perturb(
        self,
        particles: Union[List[int], int] = None,
        type: Union[List[int],int] = None,
        magnitude=10 ** -5,
        xdir=0,
        ydir=0,
        zdir=0,
    ) -> int:
        """
        Add a force to EITHER a type or list of types of particles OR a particle or list of particles in the simulation
        Pass in EITHER:
        - particles: A particle or list of particles to perturb
        OR
        - type: A type or list of types of particle to perturb (output from methods like "add_grains")

        And:
        - magnitude: magnitude of the perturbative force
        - xdir, ydir, zdir: direction of the perturbative force
        """
        if not particles is None and not type is None:
            raise Exception("You can perturb EITHER a particle OR a type of particle")
        if particles is None and type is None:
            raise Exception("You must perturb either a particle or a type of particle")
        
        if not type is None:
            if isinstance(type,int):
                group = f"type {type}"
            else: 
                group = f"type " + " ".join([str(typ) for typ in type])
        else:
            if isinstance(particles,int):
                group = f"id {particles}"
            else: 
                group = f"id " + " ".join([str(part) for part in particles])

        self._perturb_iter += 1
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\ngroup group_perturb_{self._perturb_iter} " + group + "\n")
            f.write(
                f"fix fix_perturb_{self._perturb_iter} group_perturb_{self._perturb_iter} gravity {magnitude} vector {xdir} {ydir} {zdir}\n"
            )

        return self._perturb_iter

    def add_viscosity(self, value, type: Union[List[int],int] = None, particles: Union[List[int],int] = None) -> int:
        """
        Add viscosity to all atoms or a set of types of atoms or a set of particles of atoms
        value - The strength of the viscosity
        type - The type or list of types of particles to which you want to add viscosity
        particles - The particle or list of particles to which you want to add viscosity

        If you pass in neither particles nor type, all particles in the simulation get viscous
        Returns: The number of the viscosity for if you want to remove this later
        """

        if not type is None and not particles is None:
            raise Exception("Can either add viscosity to types of particles or a list of particles, not both")
        
        if not type is None:
            if isinstance(type,int):
                group = f"type {type}"
            else: 
                group = f"type " + " ".join([str(typ) for typ in type])
        elif not particles is None:
            if isinstance(particles,int):
                group = f"id {particles}"
            else: 
                group = f"id " + " ".join([str(part) for part in particles])
        else: group = "all"

        self._visc_iter += 1
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write("\n")
            if group != "all":
                f.write(f"group group_viscosity_{self._visc_iter} " + group + "\n")
                f.write(
                    f"fix fix_viscosity_{self._visc_iter} group_viscosity_{self._visc_iter} viscous {value}\n"
                )
            else:
                f.write(
                    f"fix fix_viscosity_{self._visc_iter} all viscous {value}\n"
                )


        return self._visc_iter

    def remove_something(self, thing: str, number_of_that_thing: Union[List[int],int] = None):
        """
        This removes something from the simulation. Inputs:
        - thing: The kind of thing you want to remove, currently supported things are "viscosity", "particles",
            "type", "bonds", "angles" "gravity", "perturbation", "move", and "walls"
        - number_of_that_thing: The number of the thing you want to remove

        If thing is ___ then number_of_that_thing is ___ :
        - viscosity: the viscosity number
        - particles: a list of particle ids
        - type: a type of particles
        - bonds: if you pass in an int, number_of_that_thing is a type of bond. If you pass in a list, number_of_that_thing is a list of particles,
            the bonds between which will be deleted for example, if you pass in [1,2,3], all of the angles between particles 1,2, and 3 will be deleted.
            Be warned however, if particle 4 is bonded to particle 3, that bond will not be removed.
        - angles: if you pass in an int, number_of_that_thing is a type of angle. If you pass in a list, number_of_that_thing is a list of particles,
            the angles between which will be deleted, for example, if you pass in [1,2,3], all of the angles between particles 1,2, and 3 will be deleted.
            Be warned however, if particle 4 is angled to particle 3, that angle will not be removed.

        NOTE: If you
        """
        if thing not in [
            "viscosity",
            "particles",
            "type",
            "gravity",
            "perturbation",
            "move",
            "walls",
            "angles",
            "bonds"
        ]:
            raise Exception(
                "This is not something that I can remove, I can only remove either 'viscosity', 'particles','type', 'bonds', 'angles', 'perturbation, 'move', 'walls', or 'gravity' right now"
            )

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            if thing == "viscosity":
                f.write(f"unfix fix_viscosity_{number_of_that_thing}\n")
            if thing == "gravity":
                f.write(f"unfix grav\n")
            if thing == "perturbation":
                f.write(f"unfix fix_perturb_{number_of_that_thing}\n")
            if thing == "move":
                f.write(f"unfix fix_move_{number_of_that_thing}\n")
            if thing == "type":
                index_drop = self._particles[ self._particles['type'] == number_of_that_thing ].index
                self._particles.drop(index_drop , inplace=True)
                f.write(f"group group_delete type {number_of_that_thing}\n")
                f.write(f"delete_atoms group group_delete\n")
                f.write(f"group group_delete delete\n")
            if thing == "particles":
                self._particles.drop(index=number_of_that_thing, inplace=True)
                if isinstance(number_of_that_thing,int):
                    f.write(f"group group_delete id {number_of_that_thing}\n")
                else:
                    f.write(f"group group_delete id " + " ".join(str(part) for part in number_of_that_thing) + "\n")
                f.write(f"delete_bonds group_delete multi any\n")
                f.write(f"delete_atoms group group_delete\n")
                f.write(f"group group_delete delete\n")
            if thing == "bonds":
                if isinstance(number_of_that_thing,int):
                    f.write(f"delete_bonds all bond {number_of_that_thing}\n")
                else:
                    f.write(f"group group_delete id " + " ".join(str(part) for part in number_of_that_thing) + "\n")
                    f.write(f"delete_bonds group_delete bond 1*\n")
                    f.write(f"group group_delete delete\n")
            if thing == "angles":
                if isinstance(number_of_that_thing,int):
                    f.write(f"delete_bonds all angle {number_of_that_thing}\n")
                else:
                    f.write(f"group group_delete id " + " ".join(str(part) for part in number_of_that_thing) + "\n")
                    f.write(f"delete_bonds group_delete angle *\n")
                    f.write(f"group group_delete delete\n")
            if thing == "walls":
                if self._walls[number_of_that_thing - 1][1]:
                    f.write(f"unfix fix_walls_{number_of_that_thing}\n")
                f.write(f"region region_{number_of_that_thing} delete\n")
                if self._walls[number_of_that_thing - 1][0]:
                    f.write(f"group group_wall_{number_of_that_thing} delete\n")
        pass

    def custom(self, statement: str):
        """
        Add a custom statement to the simulation
        - statement: The line that you want to add
        """
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(statement + "\n")

        pass

    def design_thermo(
        self,
        thermo_every_this_many_seconds: float = None,
        thermo_list: List[str] = None,
    ):
        """
        Design the mid-simulation output to the log.lammps file. Inputs:
        - thermo_every_this_many_seconds: Output the thermo every this many seconds
        - thermo_list: A list of things you want in the thermo output
        """
        self._thermo_step = round(thermo_every_this_many_seconds / self._timestep)
        if not thermo_list is None:
            self._thermo_list = thermo_list
        else:
            self._thermo_list = ["step", "atoms", "time"]
        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write("\n")
            f.write("thermo_style custom " + " ".join(self._thermo_list) + "\n")
            f.write("thermo_modify lost warn\n")
            f.write(f"thermo {self._thermo_step}\n")

        pass

    def design_dump_files(
        self, dump_file_every_this_many_seconds: float, dump_list: List[str] = None
    ):
        """
        Also mention that if we have already made the dump files, this just changes how often we print them

        This command makes it such that we output dump files, and it also decides what will be in those dump files. 
        If we have already called this function once in the script, then this does not re-design the dump files,
        rather, it simply allows you to change how often those dump files will be output.
        Inputs:
        - dump_file_every_this_many_seconds: This sets how many seconds the simulation waits between printing
            dump files.
        - dump_list: This determines what will be in the dump file.
        
        The way that I have this set up, we always include the id, type, and radius of the particles, as well
        as their x y z positions. If you pass in an empty list here then that is all you'll get in your dump files.
        If you leave this set to None, we will also include the x, y, and z force components, as well as the bending
        and stretching energy. You can also additionally pass in 'contacts', which will give you the number of contacts
        per particle, 'pressure', which will output the pressure on each particle, or anything else that LAMMPS will
        accept in a custom dump file (https://docs.lammps.org/dump.html, the section on 'custom' attributes)
        """

        if self._dump_file_every_this_many_seconds:
            self._dump_file_every_this_many_seconds = dump_file_every_this_many_seconds
        
        else:
            self._dump_file_every_this_many_seconds = dump_file_every_this_many_seconds
            if not dump_list is None:
                self._dump_list = dump_list
            else:
                self._dump_list = ["fx", "fy", "fz", "bending_energy", "stretching_energy"]

            with open(os.path.join(self._path, "in.main_file"), "a") as f:
                f.write("\n")
                if "contacts" in self._dump_list:
                    f.write("compute contacts all contact/atom\n")
                    self._dump_list.remove("contacts")
                    self._dump_list.append("v_contacts")
                if "pressure" in self._dump_list:
                    f.write("compute stress all stress/atom NULL\n")
                    f.write("compute contacts all contact/atom\n")
                    f.write(
                        "variable pressure atom 2*(c_stress[1]+c_stress[2]+c_stress[3])/(c_contacts+.001)\n"
                    )
                    self._dump_list.remove("pressure")
                    self._dump_list.append("v_pressure")
                if "bending_energy" in self._dump_list:
                    f.write("compute bendingE all pe/atom angle\n")
                    self._dump_list.remove("bending_energy")
                    self._dump_list.append("c_bendingE")
                if "stretching_energy" in self._dump_list:
                    f.write("compute stretchingE all pe/atom bond\n")
                    self._dump_list.remove("stretching_energy")
                    self._dump_list.append("c_stretchingE")

                f.write(
                    f"dump pump all custom 1 out*.dump id type radius x y z "
                    + " ".join(str(item) for item in self._dump_list)
                    + "\n"
                )
                f.write("dump_modify pump pad 11\n")
        pass

    def run_simulation(self, time: float, timestep: float = None):
        """
        Run the simulation for a given amount of time. Inputs:
        - time: Amount of time to run the simulation for
        - timestep: The intended timestep of the simulation.

        The auto-generated timestep is currently in-production (TODO), when it is done, read the following:
        A note on this, a timestep is already estimated automatically via this script, and if 
        you don't pass anything into the timestep value, the simulation will use this automatically 
        calculated timestep. However, there will probably be plenty of situations where you will 
        need to change this timestep manually. If the timestep that the script automatically calculates
        is too small, the simulations will take a long time. If instead the auto-generated timestep is
        too large, the simulation will be unstable, and you might have your particles flying all over the
        place (resulting in errors where lammps tells you that it can no longer find the atoms in a bond or angle).
        If either of these things are happening to you, you might want to manually change the timestep. In that
        case, the auto-generated timestep can be a good jumping-off point!

        Another note, LAMMPS will not let you re-set the timestep if you have already run some of the simulation, 
        and then have applied a fix move. That is, if you have already simulated something -- called the 
        run_simulation() method -- and then you call the move() method, LAMMPS freaks out, because it bases its
        movement calculation on the timestep. If you change the timestep, that calculation is now out of whack.

        Based on this, if you have both:
        - Already run some of the simulation, AND
        - Applied a fix move
        Then this method will not allow you to reset the timestep
        """
        if not self._timestep:
            raise Exception("It seems like we don't have any particles to simulate!")

        if timestep is None:
            timestep = self._timestep

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            if self._have_run and self._move_iter>0:
                print("\n######")
                print("Cant reset timestep after sim has run and `fix move` has been performed (move() has been called)")
                print("Ignoring the timestep and running simulation with the last timestep that was set")
                print("######\n")
                number_of_timesteps_run = round(time / self._manual_timestep)
                number_of_timesteps_dump = round(self._dump_file_every_this_many_seconds / self._manual_timestep)
            else:
                f.write(f"\ntimestep {timestep}\n")
                number_of_timesteps_run = round(time / timestep)
                number_of_timesteps_dump = round(self._dump_file_every_this_many_seconds / timestep)
                self._manual_timestep = timestep
            if self._dump_file_every_this_many_seconds:
                f.write(
                    f"dump_modify pump every {number_of_timesteps_dump}\n"
                )
            f.write(f"run {number_of_timesteps_run}\n")
            self._have_run = True

        pass

    def manually_edit_timestep(self,timestep: float):
        """
        First thing you need to know about this method, you usually set the timestep when you run run_simulation(),
        so more often than not, this method is not needed. However, there are some rare cases where you might
        want to set this manually. For example, lets say I'm about to apply a fix move (call the move() method),
        but I also want to change the timestep before this happens. I won't be able to change the timestep after
        the fix move because of what I mention in the desription of run_simulation(). Therefore, you can use
        this function to slip in a timestep change right before the fix move.
        """

        if self._have_run and self._move_iter>0:
            raise Exception("Cant change the timestep if you have started the simulation and applied a fix move")

        with open(os.path.join(self._path, "in.main_file"), "a") as f:
            f.write(f"\ntimestep {timestep}\n")
            self._manual_timestep = timestep
    
        pass


def _dist(coords, i, j):
    return ((coords[i, 0] - coords[j, 0])**2 + (coords[i, 1] - coords[j, 1])**2 + (coords[i, 2] - coords[j, 2])**2)**0.5
