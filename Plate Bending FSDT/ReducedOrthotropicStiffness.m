function [Qbar, Pbar]  = ReducedOrthotropicStiffness(E1,E2,NU12,G12,G23,G13,theta)
% Reduced Stiffness  This function returns the reduced stiffness matrix for
% fiber reinforced materials.

NU21=NU12*E2/E1;
Q = [E1/(1-NU12*NU21) NU12*E2/(1-NU12*NU21) 0 ;NU12*E2/(1-NU12*NU21) E2/(1-NU12*NU21)  0 ; 0 0 G12];

P = [ G23 0 ; 0 G13];

m=cos(theta*pi/180);
n=sin(theta*pi/180);
T= [m*m n*n 2*m*n;n*n m*m -2*m*n;-m*n m*n m*m-n*n]; % Stress Transformation matrix
Tinv=[m*m n*n -2*m*n;n*n m*m 2*m*n;m*n -m*n m*m-n*n];
Qbar =Tinv*Q*(Tinv)';

T1 = [ m n ; -n m];

Pbar = T1*P*T1';

  