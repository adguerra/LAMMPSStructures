# LAMMPSStructures

NOTE: Can't have more than 32 groups

## Intro

The goal of this repository is to document and provide some examples for a python package which facilitates the use of LAMMPS as a simulation tool for elastic materials. This python package does not simulate things on its own, rather, it contains a class and functions that allow you to write a set of files which can then be run with LAMMPS.

Elastic structures are traditionally simulated using a finite element software, which has the benefit that it can model complex geometries. However, often these softwares have a very hard time modeling contact and other complex boundary conditions. For instance, when roots dig into soil or when cells jam on the surface of the Extra-Cellular Matrix or when birds build a nest, there is a coupling either between many elastic bodies in contact, or between granular and elastic materials. These would all be nearly impossible to simulate using a finite element software.

On the other hand, LAMMPS, a Discrete Element simulation software, is excellent at handling contact at the expense of geometrical complexity -- almost everything is treated as a sphere. In a couple of recent papers [[1](https://arxiv.org/abs/2209.05660),[2](https://arxiv.org/abs/2210.11324),[3](https://pubs.rsc.org/en/content/articlelanding/2021/sm/d1sm00787d),[4],[5]], I and others have shown that it is possible to get the best of both worlds. By "gluing" spheres together to form some simple geometries that can be endowed with elastic properties, we can have contact and elasticity in the same simulation. This repository, I hope, should contain everything you need to start doing this yourself.

## Brief Tour

We have three main directories in this repo:

 - `lammps_pypack`: This contains the python package that we will use to make the files that will allow us to simulate elastic and granular materials
 - `examples`: These are some examples of this package being used in action. Some of these examples are systems that we have studied in specific publications, but others are just examples that we thought would be intersting for future study, or which display the capacity of this tool.
 - `validation`: Here there are some simulations which we validate either theoretically or against other simulations in COMSOL. This directory is comparatively scant, just because we have also performed validation as part of some of the studies cited earlier which use this tool [[1](https://arxiv.org/abs/2209.05660),[2](https://arxiv.org/abs/2210.11324),[3](https://pubs.rsc.org/en/content/articlelanding/2021/sm/d1sm00787d),[4],[5]]

 ## How to use this tool

To begin to use this tool, clone this repository and `pip install lammpsWithPython` from within the `lammps_pypack` folder. You can then run an example file, or make a new file which sets up some simulation that you are intersted. Any time you use this tool, it will create a new folder with a name that you specify, and within that folder, it will write some files that LAMMPS will need to run a simulation. It will always write a file called `in.main_file`, which is the main input LAMMPS file that you will eventually run. This sets up the simulation, determines the dimensionality and size of your simulation box, etc. The tool may also create some secondary files which will be accessed by `in.main_file`. For example, if you want to insert a grain (a sphere) in your simulation box, you will call the `add_grains` command, and input to that command, among other things, the coordinates of those grains. When you do that, a new file will be written, called `grains_*.txt` where the asterisk is some integer (namely, 1 if this is the first time you have inserted grains, 2 if it is the second, etc.). Calling the `add_grains` command also writes a line to `in.main_file` (`include grains_*.txt`) which makes LAMMPS look in the `grains_*.txt` file for the locations of the grains. You therefore have to keep all of these files together, which is why they are all written together in a new folder every time you use the tool.

To run any of these files you will first need to download and LAMMPS. Alternatively, there are also some pre-built versions of LAMMPS. Once you have LAMMPS installed, you can go into the folder which contains `in.main_file` and input it to LAMMPS (something like `lmp_serial -i in.main_file`). This tool is currently set to run with the LAMMPS version which was released on June 23, 2022, and so if something doesn't work, let me know and I can try to adjust it! Or adjust it yourself (thats what Git is for right??).


## Some Examples

### "Cells" growing on a manifold

![cells](https://user-images.githubusercontent.com/43476955/137187691-48a16b9c-6b9f-458d-8da6-ab99c0110dac.gif)

### Wrinkles in a Developable Cone (previously undocumented)

![dcone](https://user-images.githubusercontent.com/43476955/137188365-9ddbedf4-f91a-4d99-8040-812077e63f5d.gif)

### Simulation of a knot 

![knotsgif](https://user-images.githubusercontent.com/43476955/137184828-72b2bb9d-7260-4df2-89c1-94fa9e60d076.gif)

### Frustrated growth on a cylinder

![cylinder](https://user-images.githubusercontent.com/43476955/137187100-2c82dc44-7658-4ccd-803e-c660a4c724d9.gif)

