function [Amat,Bmat,Dma] = ABD()

material =1;
materialtype = 1;

if materialtype == 1
h = b1/4;              % Thickness of plate
elseif materialtype == 2   %Isotropic
h = b1/100;
end
%% PROPERTIES OF ISOTROPIC PLATE
E = 200*10^9;
NU = 0.25;
D_iso = E*h^3/(12*(1-NU^2));

%% PROPERTIES OF LAMINATES
theta = [0 90 90 0];    % Orientation of Laminates in Laminate
NL = length(theta) ;    % Number of Layers of Laminate
if material == 1
E1(1:NL) = 75.994*10^9;     
E2(1:NL) = E1(1:NL)/25;
G12(1:NL) = 0.5*E2(1:NL);
NU12(1:NL) = 0.25;
G13(1:NL) = 0.5*E2(1:NL);
G23(1:NL) = 0.2*E2(1:NL);
elseif material == 2
E1(1:NL) = 75.994*10^9;     
E2(1:NL) = E1(1:NL)/40;
G12(1:NL) = 0.5*E2(1:NL);
NU12(1:NL) = 0.25;
G13(1:NL) = 0.5*E2(1:NL);
G23(1:NL) = 0.5*E2(1:NL);
elseif material == 3   % Graphite Epoxy
E1(1:NL) = 40*10^9;     
E2(1:NL) = 1*10^9;
G12(1:NL) = 0.6*10^9;
NU12(1:NL) = 0.25;
G13(1:NL) = 0.6*10^9;
G23(1:NL) = 0.5*10^9;
elseif material == 4     % Boron Epoxy
E1(1:NL) = 206.9*10^9;     
E2(1:NL) = E1(1:NL)/10;
G12(1:NL) = 6.9*10^9;
NU12(1:NL) = 0.3;
G13(1:NL) = 6.9*10^9;
G23(1:NL) = 4.14*10^9;
elseif material == 5
E1(1:NL) = 75.994*10^9;     
E2(1:NL) = E1(1:NL)/40;
G12(1:NL) = 0.6*E2(1:NL);
NU12(1:NL) = 0.25;
G13(1:NL) = 0.6*E2(1:NL);
G23(1:NL) = 0.5*E2(1:NL);
end


%%  ------------------   A, B and D, E, F, H and L MATRICES ------------------------------%%
Amat = zeros(3,3);
Bmat = zeros(3,3);
Dmat = zeros(3,3);
Pmat = zeros(2,2);
zp = layer_position(NL,h);    % Position of top surface of lamina w.r.t mid surface
for nl = 1:NL
    if materialtype == 1 
    [Qbar{1,nl}, Pbar{1,nl}] = ReducedOrthotropicStiffness(E1(nl),E2(nl),NU12(nl),G12(nl),G23(nl),G13(nl),theta(nl));
    elseif materialtype ==2
    [Qbar{1,nl}, Pbar{1,nl}] = IsotropicStiffness(E,NU);
    end
    
    Dmat  =  Dmat + 1/3*(zp(nl)^3 - zp(nl+1)^3)*Qbar{1,nl};
    
    Bmat  =  Bmat + 1/2*(zp(nl)^2 - zp(nl+1)^2)*Qbar{1,nl};
    
    Amat  = Amat +  (zp(nl) - zp(nl+1))*Qbar{1,nl};
    
    Pmat = Pmat + (zp(nl) - zp(nl+1))*Pbar{1,nl};
    
    
end