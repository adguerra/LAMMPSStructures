# Documentation

Here I'll split all of the functions that we have access to from `lammpsWithPython.lammps_object` into a couple of categories, go through each, and give some description and examples. The general point of having all of these commands is so that you never have to deal with LAMMPS directly, that is, you can use these python commands to cook up entire LAMMPS files and never have to think about LAMMPS as a language. There will be a few exceptions to this:

- When we execute some of these commands, specifically the `connect_particles_to_elastic`, `add_walls`, `move`, `perturb`, and `add_viscosity` commands, we will be sectioning off `group`s of particles. That is, we will only do stuff (like add a force) to a certain subset of particles. LAMMPS has a fun funky thing where it does not let us create more than 32 groups in a single simulation. This will threaten to limit you at times, for example, if you want to add 100 beams to a simulation and perturb each of them in a different direction, you can't apply a `move` to each of them, you'll have to come up with some thoughtful solutions at times (I have some in some of the example files). But hopefully, in general, this won't hold you back much.

- When you create a granular or cohesive potential between particles and grains, you will have to go to the LAMMPS documentation to understand what exactly all of the variables mean. I have those linked below

- When you create a "region" using this code, for example, if you want to restrict the motion of a bunch of grains to the inside of a cylinder, or you want to make a wall that something will run into, you'll have to go to the LAMMPS documentation to understand how to specify the geometry of those regions.

There also might be other choices that I have made that will hold you back. For example, if you feel like you wish you could pass a variable into the `move` command, that is, if you want to give some set of particles in your simulation a variable velocity, you'll be mad at me because I don't let you do that. One solution would be running the simulation iteratively, by calling `run_simulation` over and over again with small run times, each time slightly adjusting the velocity of the particles in question. You could get a lot of things to happen this way, for example, if you want to cause some elastic object to grow, you could do the same iterative running solution, but this time adjust the bond lengths every couple of nanoseconds or whatever. In this and similar cases you might also learn a little LAMMPS and use the `custom` command to make your own variable forces or velocities.

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

### `add_grains`

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

### `add_beam`

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
    
    - The ids of the particles, bonds, and angles that make up this beam
    
    Examples:
    ```python
    # Adds a vertical beam with 30 particles
    beam, _, _ = sim.add_beam(30, np.array([0.1,0,-0.05]),np.array([0.1,0,0.05]),0.005,10 ** 3,0.1)
    ```

### `add_loop`

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
    
    - The ids of the particles, bonds, and angles that make up this loop
    
    Examples:
    ```python
    # Adds a loop with 50 particles, center 000, radius radius, normal in the +z direction, etc.
    lt,_,_ = sim.add_loop(50,np.array([0,0,0]),radius, np.array([0,0,1]),h,E,density)
    ```

### `add_circular_sheet`

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
    - The ids of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a sheet with center 000, radius sheet_radius, etc.
    sheet_type, sheet_bonds, sheet_angles = sim.add_circular_sheet(np.array([0,0,0]), sheet_radius,0.004, 0.002, 10**4, 0.5)
    ```

### `add_rectangular_sheet`

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
    - The ids of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a sheet with center 000, side lengths side_length, etc.
    sheet_type, sheet_bonds, sheet_angles = sim.add_rectangular_sheet(np.array([0,0,0]), side_length, side_length,mesh_diameter, 0.002, 10 ** 4, 0.5)
    ```



### `add_cylindrical_sheet`

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
    - The ids of the particles, bonds, and angles that make up this sheet

    Examples:
    ```python
    # Adds a horizontal cylindrical sheet
    sheet_type, sheet_bonds, sheet_angles = sim.add_cylindrical_sheet(np.array([-0.075,0,0]), radius,height,np.array([1,0,0]),mesh_diameter, 0.002, 2 * 10 ** 4, 0.5)
    ```

# Make things interact

Now that you have stuff in your simulation, you may want that stuff to interact with each other or other things like walls. There are a couple of different kinds of interaction potentials that I thought might be useful to be able to deal with. Some of these I will expect to be used a lot, like `turn_on_granular` and `cohesive_potential`, `connect_particles_to_elastic`, and `add_walls`. Others, like `construct_many_bonds` and `construct_many_angles` I don't expect you to use that much. The main reason that these are in here is that they are called whenever you make some elastic geometry. For example, if you call `add_beam`, it will first place a bunch of particles in a line by calling `add_grains`, and then it will connect those grains using bonds (by calling `construct_many_bonds`) to give the beam stretching stiffness, and also it will add angles (by calling `construct_many_angles`) to give the beam bending stiffness. But I decided to include them here and make them callable themselves because who knows, maybe you'll want to manually construct bonds and angles sometime!

### `turn_on_granular_potential`

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

### `turn_on_cohesive_potential`

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


### `connect_particles_to_elastic`

- Create bonds between a set of particles and an elastic material. This will set the rest length of the bonds equal to the distance, at the time of the command, between each particle and the surface. Takes in:

    EITHER

    - `particles`: a list of particles

    OR 
    
    - `type`: a type of particles

    and:
    
    - `elastic_type`: The type of the elastic
    - `stiffness`: The stiffness of each bond
    - `n_bonds_per_particle`: How many connections you want between each particle and the elastic. This is to restrict degrees of freedom. This can be 1 2 or 3.

    Note that this cannot be undone once it is done, unless you delete the particles or the surface

    Examples:
    ```python
    # Turns on a tripod harmonic bond between a sheet of type sheet_type, and grains of type grain_type, with a stiffness in each bond of 1
    sim.connect_particles_to_elastic(sheet_type, 1, type = grains_type)
    ```


### `add_walls`

- Make a lammps "region", and optionally make the walls of that region have a granular potential with a type of particles. This can either be a simple contact potential, which repels particles which are overlapping, or it can be a super complicated potential which adds all sorts of granular mechanics. To make the region, you have to add the:

    - `region_details`: Everything that is required to define a lammps region other than the region name. For example, if I wanted to define a cylinder which was aligned with the z axis and whose center was at x=0, y=0 with a radius of 0.5 and which went from z=-0.2 to z=0.2 then I would have `region_details = cylinder z 0 0 0.5 -0.2 0.2`. The full set of options are detailed here: https://docs.lammps.org/region.html

    Optionally:

    - If you pass in a `youngs_modulus`, the walls of the region that you are defining will become stiff to some or all of the grains. Further options are:
    - If you pass in neither a `particle` or list of particle ids (particles) nor a `type` or list of types of particles,
        ALL particles in the simulation will interact with the walls
    - If you pass in either `particles`, or `type`, the particles or type of particles that you pass in will interact with the walls
    - If you want to add additional physics, such that these particles actually interact with the wall like some sort of granular material, then you should input a `hardcore_dict` which contains:
        - `restitution`: The restitution coefficient
        - `poissons`: The Poisson's ratio of the interaction
        - `xscaling`: A scalar multiplier for the tangential damping
        - `coeffric`: The coefficient of sliding friction
        - `gammar`: The degree of rolling damping
        - `rolfric`: The coefficient of rolling friction

    Outputs:

    - The id of the wall
        
    I'd check out the lammps documentation for the "pair_style granular command" if you wanna be hardcore. Note that the downside to being hardcore is that it makes the simulation take much longer

    One more thing -- if you pass in nothing for the region details, this will make the walls of the simulation station hard to the particles or type that you input, with youngs modulus and hardcore dict acting as usual. Note that this cannot currently be undone.


    Examples:
    ```python
    # Makes all particles in a simulation interact with a cylindrical wall with a youngs_modulus of E_grains
    hold2 = sim.add_walls(f"cylinder z 0 0 {column_radius} {-sim._height/2} {sim._height/2}", youngs_modulus = E_grains)
    ```
    


### `construct_many_bonds`

- Add harmonic bonds between particles. Takes in:

    - `tuplets`: An Nx2 np.array where the rows are pairs of particles that you want to make a bond between
    - `stiffness`: The stiffness of these bonds
    - `rest_length`: The rest length of the bonds

    Outputs:

    - The id of these bonds


### `construct_many_angles`

- Add cosine angles between particles. Takes in:

    - `triplets`: An Nx3 np.array where the rows are triplets of particles that you want to make an angle between
    - `stiffness`: The stiffness of these angles
    
    Outputs:
    - The id of these angles


# Do Stuff to stuff

Now that we have added grains and elastic materials to the simulation and made them interact, we might want to do something to them. Here we have four different commands, one where we can manually `move` particles, two where we add forces to particles that are operationally very similar but conceptually different (`add_gravity` and `perturb`), and `add_viscosity` which provides a drag force to the particles.


### `move`

- Move a set of particles at some velocity. Takes in 

    EITHER:

    - `particles`: A partice or list of particles to move

    OR

    - `type`: A type or list of types of particle to move (output from methods like `add_grains`)

    And:
    - `xvel`, `yvel`, `zvel`: The velocity in the x, y, and z directions to move that particle

    Note: If you pass in 'None' to either xvel, yvel, or zvel, or leave them blank, those velocities will not be mandated, and the particles will be free to move in those directions

    Examples:
    ```python
    # Move particles of type grains_type1 a velocity of 0 in the x and y direction, and -velocity in the z direction
    move1 = sim.move(type = grains_type1,xvel = 0, yvel = 0,zvel = -velocity)
    ```

### `add_gravity`

- Add gravity to the simulation. Takes in:

    Optionally:
    - `magnitude`: magnitude of the gravitational force
    - `xdir`, `ydir`, `zdir`: direction of the gravitational force

    The default values are:
    - `magnitude`: 9.8
    - `xdir`: 0
    - `ydir`: 0
    - `zdir`: -1

    Examples:
    ```python
    # Add a regular 'ol gravitational force to all particles in the simulation
    sim.add_gravity()
    ```

### `perturb`

- Add a force to particles in the simulation. Takes in:
    
    EITHER:
    - `particles`: A particle or list of particles to perturb

    OR
    - `type`: A type or list of types of particle to perturb (output from methods like "add_grains")

    And:
    - `magnitude`: magnitude of the perturbative force
    - `xdir`, `ydir`, `zdir`: direction of the perturbative force

    Default values are:
    - `magnitude`: 10 ** -5
    - `xdir`: 0
    - `ydir`: 0
    - `zdir`: 0


### `add_viscosity`

- Add viscosity to all atoms or a set of types of atoms or a set of particles of atoms. Takes in:
    - `value`: The strength of the viscosity

    And EITHER
    - `type`: The type or list of types of particles to which you want to add viscosity

    OR
    - `particles`: The particle or list of particles to which you want to add viscosity

    If you pass in neither particles nor type, all particles in the simulation get viscous

    Outputs:
    - The id of the viscosity for if you want to remove this later

    Examples:
    ```python
    # Add a viscosity to all particles in the simulation with coefficient "viscosity"
    sim.add_viscosity(viscosity)
    ```

# Output

There are two main ways that you can output information from LAMMPS during the simulation run. One is whats called a `thermo` which is short for thermodynamical properties. At every timestep, some amount of information (which you choose) can be outputted to the `log.lammps` file as well as the console. I don't use this very much, but you could for example use it to output the number of particles at every 1000 timesteps to make sure no particles have escaped while the sim is running. Another output which I use all the time for every simulation are `dump` files. They can be input into visualization tools like OVITO to let you see what went on in your simulation. NOTE that you can only use OVITO to see particles, you can't see, for example, walls that you add in using the `add_walls` command.

### `design_thermo`

- Design the mid-simulation output to the log.lammps file. Takes in:
    - `thermo_every_this_many_seconds`: Output the thermo every this many seconds
    - `thermo_list`: A list of things you want in the thermo output

### `design_dump_files`

- This command makes it such that we output dump files, and it also decides what will be in those dump files. If we have already called this function once in the script, then this does not re-design the dump files, rather, it simply allows you to change how often those dump files will be output. Takes in:
    - `dump_file_every_this_many_seconds`: This sets how many seconds the simulation waits between printing dump files.
    - `dump_list`: This determines what will be in the dump file.
    
    The way that I have this set up, we always include the id, type, and radius of the particles, as well as their x y z positions. If you pass in an empty list here then that is all you'll get in your dump files. If you leave this set to None, we will also include the x, y, and z force components, as well as the bending and stretching energy. You can also additionally pass in 'contacts', which will give you the number of contacts per particle, 'pressure', which will output the pressure on each particle, or anything else that LAMMPS will accept in a custom dump file (https://docs.lammps.org/dump.html, the section on 'custom' attributes). 

    Examples:
    ```python
    # Output the default dump files every 0.01 seconds
    sim.design_dump_files(0.01)
    ```

# Misc

### `run_simulation`

- Run the simulation for a given amount of time. Takes in:
    - `time`: Amount of time to run the simulation for
    - `timestep`: The intended timestep of the simulation.

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
    - Already run some of the simulation, 
    
    AND
    - Applied a fix move
    Then this method will not allow you to reset the timestep

### `custom`

- Add a custom statement to the simulation. Takes in:
    - `statement`: The line that you want to add

### `remove_something`

- This removes something from the simulation. Takes in:
    - `thing`: The kind of thing you want to remove, currently supported things are `viscosity`, `particles`,
        `type`, `bonds`, `angles` `gravity`, `perturbation`, `move`, and `walls`
    - `id_of_that_thing`: The id of the thing you want to remove

    If `thing` is _____ then `id_of_that_thing` is ____ :
    - `viscosity`: the viscosity id
    - `particles`: a list of particle ids
    - `type`: a type of particles
    - `bonds`: if you pass in an int, id_of_that_thing is a type of bond. If you pass in a list, id_of_that_thing is a list of particles,
        the bonds between which will be deleted for example, if you pass in [1,2,3], all of the angles between particles 1,2, and 3 will be deleted.
        Be warned however, if particle 4 is bonded to particle 3, that bond will not be removed.
    - `angles`: if you pass in an int, id_of_that_thing is a type of angle. If you pass in a list, id_of_that_thing is a list of particles,
        the angles between which will be deleted, for example, if you pass in [1,2,3], all of the angles between particles 1,2, and 3 will be deleted.
        Be warned however, if particle 4 is angled to particle 3, that angle will not be removed.

    Examples:
    ```python
    # Removes the gravity from the simulation
    sim.remove_something("gravity")
    # Removes a wall with id `wall1` from the simulation
    sim.remove_something("walls",wall1)
    ```

### `manually_edit_timestep`

- First thing you need to know about this method, you usually set the timestep when you run run_simulation(), so more often than not, this method is not needed. However, there are some rare cases where you might want to set this manually. For example, lets say I'm about to apply a fix move (call the `move`), but I also want to change the timestep before this happens. I won't be able to change the timestep after the fix move because of what I mention in the desription of run_simulation(). Therefore, you can use this function to slip in a timestep change right before the fix move. Takes in:
    - `timestep`: the timestep that you want to set

