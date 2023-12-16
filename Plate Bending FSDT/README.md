# Bending Analysis of Fiber-Reinforced Composite
### Main file - platefsdtfem.m
Takes plate dimensions, material properties, loading conditions, and mesh parameters ( number of nodes, degree of freedom per node ) to give equilibirium deformation of the plate under loading

## Functions
### gausspoints()
Returns the weights and points for numerical integration using the Gauss quadrature algorithm. Takes the order of integration to be performed ( can go up to order 4 ( i.e. 4x4 points )) as input.

### planestressassemblyr()
Constructs the global stiffness matrix for the loaded plate using the element-level stiffness matrix
