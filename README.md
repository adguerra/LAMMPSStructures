# LAMMPSStructures

The goal of this repository is to provide some examples and source code for the use of LAMMPS as a simulation tool for elastic materials. These materials are traditionally simulated using a finite element software, which has the benefit that it can model complex geometries, but often these softwares have a very hard time modeling contact and other complex boundary conditions. For instance, when roots dig into soil or when cells jam on the surface of the Extra-Cellular Matrix or when birds build a nest, there is a coupling either between many elastic bodies in contact, or between granular and elastic materials. These would all be nearly impossible to simulate using a finite element software.

On the other hand, LAMMPS, a Discrete Element simulation software, is excellent at handling contact at the expense of geometrical complexity -- almost everything is treated as a sphere. Here we show how we have used the tools in LAMMPS to "glue" spheres into elastic bodies to get the best of both worlds. We give some examples for how we have used this tool so far, and provide some example code to look through. Eventually we would like to turn this repository into a tool, either a general LAMMPS input file which is easy to use, or some sort of GUI or MATLAB script where one could input geometry and load and then run a simulation of their elastic structure

## Brief Tour

In this main folder we have a first draft of a PDF of the derivation of all of the interatomic forces -- PairPotentialDerivations.pdf. Until further notice we recommend instead referring to [this publication](https://pubs.rsc.org/en/content/articlelanding/2021/sm/d1sm00787d) for a more thorough, well-edited, and specific derivation of the pair potentials for 1-d elastic structures.

We split our examples into two categories:

 - 1d structures, which are "grains" or "particles" arranged in a line to form an elastic beam
 - 2d structures, which are particles arranged in a triangular lattice to form an elastic sheet

In each of these directories I will have some verification (probably in COMSOL, or with experiments from our lab) for some different kind of structures, videos of some of our simulations and some heavily commented example files. 

The most heavily commented, and also the simplest, LAMMPS input file is LAMMPSStructures/1d/ValidateCylindricalBeams/in.cylindricalValidation, and the other 1d files are commented where they differ from that file. If you are new to this repository or to LAMMPS I would start there. The most heavily commented 2d file is LAMMPSStructures/2d/REPLACEWHENYOUVALIDATE2D.

 To run any of these files you will need to download LAMMPS

List and description of files: 

1d
 -  ValidateCylindricalBeams: Code to perform a test of the buckling strain of cylindrical beams as well as comparrison to COMSOL
 -  ValidateRectangularBeams: Same as above but with rectangular beams
 -  ValidateLoops: Code to test the shape of compressed rings and comparrison to theory
 -  2BeamCompetitionRectangular: Code which compresses two beams next to one another
 -  ColumnWithLoops: Code for a column of grains confined with elastic loops

2d
 - Still working... 


## Some Examples

### "Cells" growing on a manifold

![cells](https://user-images.githubusercontent.com/43476955/137187691-48a16b9c-6b9f-458d-8da6-ab99c0110dac.gif)

### Wrinkles in a Developable Cone (previously undocumented)

![dcone](https://user-images.githubusercontent.com/43476955/137188365-9ddbedf4-f91a-4d99-8040-812077e63f5d.gif)

### Simulation of a knot 

![knotsgif](https://user-images.githubusercontent.com/43476955/137184828-72b2bb9d-7260-4df2-89c1-94fa9e60d076.gif)

### Frustrated growth on a cylinder

![cylinder](https://user-images.githubusercontent.com/43476955/137187100-2c82dc44-7658-4ccd-803e-c660a4c724d9.gif)

