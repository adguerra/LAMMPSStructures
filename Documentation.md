NOTE: Can't have more than 32 groups


# Documentation

Here I'll split all of the functions that we have access to from `lammpsWithPython.lammps_object` into a couple of categories, go through each, and give some description and examples.

# Startup

When you first start up a simulation, you will create an instance of the `Simulation` class by doing something like:
```python
# Example 1
sim = Simulation("simulation_dir", 3, 0.5, 0.5, 0.2)
# Example 2
sim = Simulation("simulation_dir", 3, 0.5, 0.5, 0.2, y_bound = "s", z_bound = "p")
```
This specific example opens a new directory called `simulation_dir` and starts to write a file `in.main_file` within that folder. This `in.main_file` will be the LAMMPS file that you will eventually run to run your simulation. Further functions that you call on the `sim` object that you created here will add to this file, and write other files that will be referenced by `in.main_file`. This initial command is also used to set some features of the simulation that you are about to create, for example, in Example 1 above, we have made a 3 dimensional simulation, with a width (x-dim) of 0.5, a length (y-dim) of 0.5, and a height (z-dim) of 0.2 with fixed boundary conditions (the simulation box is always centered at [0,0,0]). In Example 2, we have done the same thing except the y-boundaries of the box are shrink-wrapped, and the z-boundaries are periodic (see the LAMMPS documentation for further descriptions of these bcs). All of these features are written to `in.main_file` along with a bunch of other set-up stuff that I will not describe here, but which I have shown and commented in `lammpsWithPython/lammps_object.py`.

# Add stuff into a simulation

Whenever you add grains or some elastic body into a simulation, this tool writes a new file, for example, if you haven't added any beams to your simulation yet and you call the `add_beams` function, it will create a file `beam_1.txt` inside of your simulation folder (the next time you add a beam the file will be called `beam_2.txt`, etc.). It also writes a line to `in.main_file` to call this other file, something like `include beam_1.txt`

## `add_grains`

- Add a bunch of individual grains to a simulation. Takes in:

    - `coords`: an Nx3 numpy array where N is the number of grains that we are trying to insert and
        the three columns are the x, y, and z coordinates of those grains
    - `diameter`: the diameter of the grains that you will be inserting
    - `density`: the density of the grains that you will be inserting
    
    Optionally:

    - `filename`: the name of the file that you want to put the create lines in. Its recommended to leave this empty, this is mostly for use by the other functions
    
    Outputs:
    
    - The type id of the particles that you are adding. This is for input into methods like "add_pair_potential" or "remove_something"
    
    Examples:
    ```python
    # This designs an Nx3 array of coordinates of grains we want to add to the simulation
    grains = np.column_stack(((np.random.rand(number_of_grains, 2) - 0.5),np.zeros(number_of_grains)))
    # We add them in at those coordinates, with some diameter, and a density of 0.01
    grains_type = sim.add_grains(grains, grain_diameter, 0.01)
    ```

## `add_beam`

- Create a beam out of particles, replete with bonds and angles to give the beam mechanical properties. Takes in:

    - `n_particles`: The number of particles that you want to make up your beam
        the three columns are the x, y, and z coordinates of those grains
    - `pos1`: A 1x3 np.array with the x, y, z coordinates of one end of the beam
    - `pos2`: A 1x3 np.array with the x, y, z coordinates of the other end of the beam
    - `geometric_thickness`: The diameter of the particles which make up the beam
    - `youngs_modulus`: The Young's modulus of the beam
    - `density`: The density of the particles which make up the beam. NOTE that this is different from
        the density of the beam itself
    
    Optionally:

    - `energetic_thickness_multiplier`: This is to change the thickness term in the stretching and bending
        energies, which changes the stretching energy linearly and the bending energy cubically
    - `stretching_multiplier`: In case you want to manually edit the stretching energy
    - bending_multiplier: In case you want to manually edit the bending energy
    
    Outputs:
    
    - The number of the particles, bonds, and angles that make up this beam
    
    Examples:
    ```python
    # Adds a vertical beam with 30 particles
    beam, _, _ = sim.add_beam(30, np.array([0.1,0,-0.05]),np.array([0.1,0,0.05]),0.005,10 ** 3,0.1)
    ```

## `add_loop`

- Create a loop out of particles, replete with bonds and angles to give the loop mechanical properties. Takes in:

    - `n_particles`: The number of particles that you want to make up your beam
    - `center`: A 1x3 np.array with the x, y, z coordinates of the center of the loop
    - `radius`: The radius of the loop
    - `normal`: A 1x3 np.array which gives the vector normal to the plane which the loop occupies,
        for example, a ring which is laying flat on a table would have a normal \porm [0, 0, 1]
    - `geometric_thickness`: The diameter of the particles which make up the beam
    - `youngs_modulus`: The Young's modulus of the beam
    - `density`: The density of the particles which make up the beam. NOTE that this is different from
        the density of the beam itself
    
    Optionally: 
    - `energetic_thickness_multiplier`: This is to change the thickness term in the stretching and bending
        energies, which changes the stretching energy linearly and the bending energy cubically
    - `stretching_multiplier`: In case you want to manually edit the stretching energy
    - `bending_multiplier`: In case you want to manually edit the bending energy
    
    Outputs:
    
    - The number of the particles, bonds, and angles that make up this loop
    
    Examples:
    ```python
    # Adds a loop with 50 particles, center 000, radius radius, normal in the +z direction, etc.
    lt,_,_ = sim.add_loop(50,np.array([0,0,0]),radius, np.array([0,0,1]),h,E,density)
    ```

## `add_circular_sheet`

- Create a circular sheet out of particles, replete with bonds and angles to give the sheet mechanical properties. For now, all sheets have a normal of [0,0,1]. I might change this later. Takes in:

    - `center`: A 1x3 np.array with the x, y, z coordinates of the center of the sheet
    - `radius`: The radius of the sheet
    - `mesh_particle_spacing`: The resting distance between mesh particles in the sheet
    - `energetic_thickness`: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
        NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
        Explicitly, for you elasticity mathematicians out there, this is h.
    - `youngs_modulus`: The Young's modulus of the sheet material
    - `density`: The density of the particles which make up the sheet. NOTE that this is different from
        the density of the sheet itself
    
    Optionally:
    - `mesh_particle_diameter`: This term changes the diameter of the mesh particles. If this is not set (left set to None),
        the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
        will look like they're right next to each other
    - `stretching_multiplier`: In case you want to manually edit the stretching energy
    - `bending_multiplier`: In case you want to manually edit the bending energy

    Outputs:
    - The number of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a sheet with center 000, radius sheet_radius, etc.
    sheet_type, sheet_bonds, sheet_angles = sim.add_circular_sheet(np.array([0,0,0]), sheet_radius,0.004, 0.002, 10**4, 0.5)
    ```

## `add_rectangular_sheet`

- Create a rectangular sheet out of particles, replete with bonds and angles to give the sheet mechanical properties. For now, all sheets have a normal of [0,0,1]. I might change this later. Takes in:

    - `center`: A 1x3 np.array with the x, y, z coordinates of the center of the sheet
    - `side_length1`: The x-dimension length of the sheet
    - `side_length2`: The y-dimension length of the sheet
    - `mesh_particle_spacing`: The resting distance between mesh particles in the sheet
    - `energetic_thickness`: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
        NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
        Explicitly, for you elasticity mathematicians out there, this is h.
    - `youngs_modulus`: The Young's modulus of the sheet material
    - `density`: The density of the particles which make up the sheet. NOTE that this is different from
        the density of the sheet itself
    
    Optionally:
    - `mesh_particle_diameter`: This term changes the diameter of the mesh particles. If this is not set (left set to None),
        the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
        will look like they're right next to each other
    - `stretching_multiplier`: In case you want to manually edit the stretching energy
    - `bending_multiplier`: In case you want to manually edit the bending energy

    Outputs:
    - The number of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a sheet with center 000, side lengths side_length, etc.
    sheet_type, sheet_bonds, sheet_angles = sim.add_rectangular_sheet(np.array([0,0,0]), side_length, side_length,mesh_diameter, 0.002, 10 ** 4, 0.5)
    ```



## `add_cylindrical_sheet`

- Create a cylindrical sheet out of particles, replete with bonds and angles to give the sheet mechanical properties. Takes in:

    - `center`: A 1x3 np.array with the x, y, z coordinates of the center of the bottom ring of the cylinder
    - `radius`: The radius of the sheet
    - `normal`: The vector which points down the center of the cylinder
    - `mesh_particle_spacing`: The resting distance between mesh particles in the sheet
    - `energetic_thickness`: This term is what is used to calculate the bending and stretching modulus of the sheet (along with the youngs modulus). 
        NOTE: changing this does not change any physical sizes in the simulation, only the energy in the bonds and angles between particles.
        Explicitly, for you elasticity mathematicians out there, this is h.
    - `youngs_modulus`: The Young's modulus of the sheet material
    - `density`: The density of the particles which make up the sheet. NOTE that this is different from
        the density of the sheet itself

    Optionally:
    - `mesh_particle_diameter`: This term changes the diameter of the mesh particles. If this is not set (left set to None),
        the mesh particle diameter will simply be set to the mesh particle spacing, such that all of the particles in the sheet
        will look like they're right next to each other
    - `stretching_multiplier`: In case you want to manually edit the stretching energy
    - `bending_multiplier`: In case you want to manually edit the bending energy

    Outputs:
    - The number of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a horizontal cylindrical sheet
    sheet_type, sheet_bonds, sheet_angles = sim.add_cylindrical_sheet(np.array([-0.075,0,0]), radius,height,np.array([1,0,0]),mesh_diameter, 0.002, 2 * 10 ** 4, 0.5)
    ```

# Make things interact

Now that you have stuff in your simulation, you may want that stuff to interact with each other or other things like walls. There are a couple of different kinds of interaction potentials that I thought might be useful to be able to deal with. Some of these I will expect to be used a lot, like `turn_on_granular` and `cohesive_potential`, `connect_particles_to_elastic`, and `add_walls`. Others, like `construct_many_bonds` and `construct_many_angles` I don't expect you to use that much. The main reason that these are in here is that they are called whenever you make some elastic geometry. For example, if you call `add_beam`, it will first place a bunch of particles in a line by calling `add_grains`, and then it will connect those grains using bonds (by calling `construct_many_bonds`) to give the beam stretching stiffness, and also it will add angles (by calling `construct_many_angles`) to give the beam bending stiffness. But I decided to include them here and make them callable themselves because who knows, maybe you'll want to manually construct bonds and angles sometime!

## `turn_on_granular_potential`

- Make two types of particles interact in a granular way. This can either be a simple contact potential, which repels particles which are overlapping, or it can be a super complicated potential which adds all sorts of granular mechanics. If you want the chill potentiall, all you need to input is:
    - `type1`, `type2`: The types of particles you want to add an interaction potential to. This is output from methods like "add_grains," "add_beam," etc. 
    - `youngs_modulus`: The youngs modulus of the particles -- this will determine how much they can overlap given some confining force

    If you don't pass in type2, it will turn on the granular potential between type1 and all other particles. If you pass in neither type1 nor type2 it will turn on the granular potential between all particles.
        
    If you want to add additional physics, such that these particles actually behave like some sort of granular material, then you should:
    - input a `hardcore_dict` which contains:
        - `restitution`: The restitution coefficient
        - `poissons`: The Poisson's ratio of the interaction
        - `xscaling`: A scalar multiplier for the tangential damping
        - `coeffric`: The coefficient of sliding friction
        - `gammar`: The degree of rolling damping
        - `rolfric`: The coefficient of rolling friction
    I'd check out the lammps documentation for the "pair_style granular command" if you wanna be hardcore. Note that the downside to being hardcore is that it makes the simulation take much longer

    Examples:
    ```python
    # Turns on a granular potential between particles of type `grains_type` and themselves, with Young's modulus E_grains, and further contact properties as shown
    contact_properties = {
        "restitution" : 0.4,
        "poissons": 0.35,
        "xscaling": 1,
        "coeffric": 0.5,
        "gammar": 1,
        "rolfric": 0.5,
    }
    sim.turn_on_granular_potential(grains_type, grains_type, E_grains, contact_properties)

    # Turns on a simple contact potential, that repels particles that are overlapping, between particles of type grains_type, and particles of type loop_type, and Youngs modulus E_grains
    sim.turn_on_granular_potential(grains_type, loop_type,E_grains)
    ```

## `turn_on_cohesive_potential`

- Make two types of particles have a cohesive (or repellant) lj potential. This takes in:
    - `type1`, `type2`: The types of particles you want to add an interaction potential to. This is output from methods like "add_grains," "add_beam," etc. 
    - `cohesive_strength`: This is epsilon in https://docs.lammps.org/pair_lj.html
    - `rest_length`: The rest length of the bonds. This is NOT sigma in https://docs.lammps.org/pair_lj.html
    - `cutoff`: If two particles get farther than "cutoff" apart, they will not be coherent any more. This can just be really big,
        or like 2.5 times rest_length, where the potential would pretty much disappear anyways

    If you don't pass in type2, it will turn on the potential between type1 and all other particles. If you pass in neither type1 nor type2 it will turn on the potential between all particles.

    Examples:
    ```python
    # Turns on a cohesive potential between particles of type grains_type and each-other
    sim.turn_on_cohesive_potential(grains_type, grains_type, grain_grain_cohesion, grain_grain_rest_length, grain_grain_rest_length * 2.5)
    ```


## `connect_particles_to_elastic`
## `add_walls`


## `construct_many_bonds`

- Add harmonic bonds between particles. Takes in

## `construct_many_angles`








# Do Stuff to stuff

## `move`
## `add_gravity`
## `perturb`
## `add_viscosity`

# Output

## `design_thermo`
## `design_dump_files`

# Misc

## `run_simulation`
## `custom`
## `remove_something`
## `manually_edit_timestep`


