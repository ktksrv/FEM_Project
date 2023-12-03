%% Clearing the Workspace
clc
clear all
close all

%% Geometrical Paramters
E = 30e6;
L = 100;
h = L/100;
wb = L/100;
Area = wb*h;
I = wb*h^3/12;


NE = 40;
NNE = 2;
NN = NNE*NE - NE + 1;
NDOF = 3;
a = L/NE;
CONN = [1:NN-1;2:NN]'; 
         

qu0 = 0;




A = [1 0 0 0; 0 1 0 0; 1 a a^2 a^3; 0 1 2*a 3*a^2];

Ainv = inv(A);
ngk = 3;
ngl = 3;
[wk,gpk] = gausspoints(ngk);
[wl,gpl] = gausspoints(ngl);

for NR = 1:2;
tol = 1e-4;
norm = 1;
qw0 = (NR-1)*1;
SKE = NNE*NDOF; % Size of Element Stiffness Matrix
SKG = NN*NDOF; % Size of Global Stiffness Matrix
d =0*ones(SKG,1);
while(norm>tol)
%% Assembly of Stiffness Matrix and Force Vector
KG = zeros(SKG);
KGT = zeros(SKG);
FG = zeros(SKG,1);
DOFM = zeros(NE,NDOF*NNE);
for i = 1:NE
KE = zeros(SKE,SKE);
KTE = zeros(SKE,SKE);
FE = zeros(SKE,1);   
    for j = 1:NNE
        nodenum = CONN(i,j);
        for k = 1:NDOF
        DOFM(i,NDOF*(j-1) + k) = NDOF*(nodenum - 1) + k;
        end
    end


for r = 1:ngk
    weight = wk(r);
    xgk = a/2*(1 + gpk(r));
   
BM = [-1/a 0 0 1/a 0 0];

BB1 = [0 0 2 6*xgk]*Ainv;

BB = [0 BB1(1) BB1(2) 0 BB1(3) BB1(4)];



% LINEAR STIFFNESS MATRIX
KEL = BM'*E*Area*BM + BB'*E*I*BB;



% NON-LINEAR STIFFNESS MATRIX

PSINL = [0 1 2*xgk 3*xgk^2]*Ainv;


G =[0 PSINL(1) PSINL(2) 0 PSINL(3) PSINL(4)];
 
dnl = d(DOFM(i,:),1);

BNL = dnl'*(G'*G);

KENL = BM'*E*Area*BNL/2 + BNL'*E*Area*BM + BNL'*E*Area*BNL/2;


% ELEMENTAL STIFFNESS MATRIX

KE1 = KEL + KENL;

%TANGENT STIFNESS MATRIX

NX = E*h*(BM + BNL/2)*dnl;

KENLT = BM'*E*Area*BNL + BNL'*E*Area*BM + BNL'*E*Area*BNL;

KSE = G'*NX*G;

KTE1 = KEL + KENLT + KSE;

     jac = a/2;
    KE = KE + weight*jac*KE1;
    KTE = KTE + weight*jac*KTE1;
end


% ELEMENTAL FORCE MATRIX 
    for s = 1:ngl
    weight = wl(s);
    xgl = a/2*(1 + gpl(s));  
PSIW1 = [1 xgl xgl^2 xgl^3]*Ainv;

PSIW = [0 PSIW1(1) PSIW1(2) 0 PSIW1(3) PSIW1(4)];

PSIU = [1-xgl/a 0 0 xgl/a 0 0];
FWE = PSIW'*qw0;

FUE = PSIU'*qu0;

FE1 = FWE + FUE;
    FE = FE + weight*jac*FE1;
    end
    [KG,KGT,FG] = beamassembly(KE,KTE,FE,NNE,NDOF,DOFM(i,:),KG,KGT,FG);
end




%% Application of Boundary Conditions
KGS = KG;
FGS = FG;
cdof = [1:2, (SKG-2:SKG-1)];
for ii = cdof
KG(ii,:) = 0;
KG(:,ii) = 0;
KG(ii,ii) = 1;
KGT(ii,:) = 0;
KGT(:,ii) = 0;
KGT(ii,ii) = 1;
FG(ii,1) = 0;
end

deltaFG = FG - KG*d ;
delta = linsolve(KGT,deltaFG);
d = d + delta;
deltamod = abs(delta);
norm = sqrt(deltamod'*deltamod)/sqrt(d'*d);
end
load(NR) = qw0;
dispw = d(2:3:end);
maxdisp(NR) = max(dispw);
wexact_linear(NR) = 5*qw0*L*L*L*L/(384*E*I);
figure(1); plot(dispw); xlabel('Node Number')
ylabel('Displacement (mm)')
hold on
end
maxdisp
wexact_linear
figure(2);plot(load,maxdisp,'*')
hold on
plot(load,wexact_linear,'o')