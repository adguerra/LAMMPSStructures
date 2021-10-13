# LAMMPSStructures

This repository is intended to showcase some 


https://user-images.githubusercontent.com/43476955/137153475-48e0dc59-a64e-490d-a1e9-bc2e374901e6.mov


## Map

In this repository I have a couple of different examples of how we make and use elastic structures in LAMMPS. I will split this up into two categories of structures:

 - 1d structures, which are "grains" or "particles" arranged in a line to form an elastic beam
 - 2d structures, which are particles arranged in a triangular lattice to form an elastic sheet

In each of these directories I will have some verification (probably in COMSOL, or with experiments from our lab) for some different kind of structures, videos of some of our simulations and some heavily commented example files. In this main folder I also have a PDF of the derivation of all of the interatomic forces -- PairPotentialDerivations.pdf.

The most heavily commented, and also the simplest, LAMMPS input file is LAMMPSStructures/1d/ValidateCylindricalBeams/in.cylindricalValidation, and the other 1d files are commented where they differ from that file. If you are new to this repository or to LAMMPS I would start there. The most heavily commented 2d file is LAMMPSStructures/2d/REPLACEWHENYOUVALIDATE2D.

 To run any of these files you will need to download LAMMPS. Hopefully there is enough information here that you could take the code and modify it to simulate your own situation :) 



List and description of files: 

1d
 -  ValidateCylindricalBeams: Code to perform a test of the buckling strain of cylindrical beams as well as comparrison to COMSOL
 -  ValidateRectangularBeams: Same as above but with rectangular beams
 -  ValidateLoops: Code to test the shape of compressed rings and comparrison to theory
 -  2BeamCompetitionRectangular: Code which compresses two beams next to one another
 -  ColumnWithLoops: Code for a column of grains confined with elastic loops

2d
 - Still working... 
