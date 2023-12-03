function [KE, FE] = beamelement(E,I,L,q0)



KE = E*I/L^3*[12   6*L   -12  6*L;...
               6*L  4*L^2 -6*L 2*L^2;...
              -12  -6*L    12  -6*L;... % Stiffness matrix in local coordinates
               6*L  2*L^2 -6*L 4*L^2];
FE = q0*L/12*[6;L;6;-L];