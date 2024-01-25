# CFD_with_FEniCS
 
 Hello ! Welcome to my CFD repository.
 
The objective of this repository is to embark on a computational odyssey where we harness the simplicity of Fenics to navigate the complexities of computational fluid dynamics (CFD).
Fenics stands out as a cutting-edge finite element library, providing a robust framework for solving partial differential equations (PDEs). Its flexibility, efficiency, and user-friendly interface makes it an ideal choice for CFD simulations. By leveraging Fenics, we aim to streamline the process of developing, implementing, and solving CFD problems, empowering both beginners and experts in the field.

I welcome contributions from the community. Whether it's bug fixes, new features, or documentation improvements, your input is valuable.

This repository is set up by folders. The title of each folder gives a brief description of what CFD model is to be solved. In general, each folder will contain 4 files. One main .py file, two .xml files containing the geometrical and mesh information and a Results folder. The Results folder contains the solution of the variables that is solved in each model and stored in a .pvd and .vtu file later to be analyzed in ParaView.
In the repository there's a folder called ZMesh that contains the code for the mesh generation through Gmsh.

It is worth noting that all the software used is open-source (Python, Fenics, Paraview & Gmsh).

FEniCS:
https://fenicsproject.org/

ParaView:
https://www.paraview.org/

Gmsh:
https://gmsh.info/
