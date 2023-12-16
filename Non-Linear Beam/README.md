# Geometrical Non-Linear Analysis of Beam using FEM
### Main file - nonlinear_beam.m 
Takes Elastic Modulus, Dimensions of the beam, No of elements or nodes, degree of freedom per node, load distribution as input and gives equilibrium state of the non-linear beam under the applied load as output

## Functions 
### gausspoints()
Returns the weights and points for numerical integration using the Gauss quadrature algorithm. Takes the order of integration to be performed ( can go up to order 4 ( i.e. 4x4 points ) as input.

### beamelement()
Returns the element level stiffness and force matrix for a beam element, given the elastic modulus, second moment of inertia, and load intensity as input.

### beamassembly() 
Constructs the global stiffness matrix for the loaded beam using the element-level stiffness matrix 


