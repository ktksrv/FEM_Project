function [Q, P] = IsotropicStiffness(E,NU)
% Reduced Isotropic Stiffness:  This function returns the reduced isotropic
% stiffness matrix for given E and NU
Q =[E/(1-(NU*NU)) NU*E/(1-NU*NU) 0 ; NU*E/(1-NU*NU) E/(1-NU*NU) 0 ;...
                                                           0 0 E/(2*(1+NU))];
G23 = E/(2*(1+NU));
G13 = G23;
P = [G23 0 ; 0 G13];                                                     
                                                       