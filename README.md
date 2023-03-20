# Minimum Cycle Generator
This code extracts connected components and minimal cycles for Homology H1 generators from spatially embedded networks. The implementation is based on the [Breadth First Search(BFS)](https://en.wikipedia.org/wiki/Breadth-first_search) algorithm and uses the [Gudhi](https://gudhi.inria.fr/) library to compute the Vietoris rips complex.
The Vietoris Rips complex is then used to compute persistent features H1 generator.
The input to the code is a spatially embedded network (network along with spatial location of nodes), and the output is the .VTP files so one can use [Paraview](https://www.paraview.org/) to visualize the results.

## Requirments
In addition to the other dependencies in requirements.txt, the implementation uses [VTK](https://vtk.org/)(version: 9) for only visualization perpose. 
The instruction to build VTK is available [here]https://vtk.org/Wiki/VTK/Configure_and_Build

