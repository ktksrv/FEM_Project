%% Plate Bending Using Reissner Mindlin Plate Theory

clc 
clear all
close all
format long

%% Plate Dimensions
a = 1;
b = 1;
h = a/100;
q0 = 1000;
%% Material Properties
Y = 210e9;
NU = 0.25;



%% Mesh Generation
gtrapx = 16;
gtrapy = 16;
numnod = (gtrapx+1)*(gtrapy+1);

NN = numnod;
NDOF = 5; % Number of degrees of freedom per node
NNE = 4;  % Number of nodes per element
NE = gtrapx*gtrapy; % Number of elements


xgtrap = a/gtrapx;   % Width of each divsion of 'a'
ygtrap = b/gtrapy;   % Height of each divsion of 'b'


% Coordinates of nodes for background mesh
[x, y] = meshgrid(0:xgtrap:a,0:ygtrap:b);
xg(1,:) = reshape(x',1,numnod);
xg(2,:) = reshape(y',1,numnod);
% Draw vertical lines
hold on
for k = 1:gtrapx+1
    nodeb = (k-1)*(gtrapy+1) + 1;
    nodet =  k*(gtrapy+1);
    xnodeb = xg(1,nodeb);  % Coordinates of bottom node
    ynodeb = xg(2,nodeb);
    xnodet = xg(1,nodet);  % Coordinates of top node
    ynodet = xg(2,nodet);
figure(1); plot([xnodeb xnodet],[ynodeb,ynodet],'k')    
end

% Draw horizontal lines
for l = 1:gtrapy+1
    nodel = l;
    noder = gtrapx*(gtrapy+1) + l ;
    xnodel = xg(1,nodel);  % Coordinates of left node
    ynodel = xg(2,nodel);
    xnoder = xg(1,noder);  % Coordinates of right node
    ynoder = xg(2,noder);
   plot([xnodel xnoder],[ynodel,ynoder],'k')
end


% SET UP CONNECTIVITY ARRAY
for j = 1:gtrapx
   for i = 1:gtrapy
   elemn = (j-1)*gtrapy + i;
   noder(elemn,1) = elemn + (j-1);
   noder(elemn,2) = noder(elemn,1) + 1;
   noder(elemn,3) = noder(elemn,1) + gtrapx + 2;
   noder(elemn,4) = noder(elemn,3)  - 1;
   
   sumx = 0;
   sumy = 0;
   for av = 1:4
       sumx = sumx + xg(1,noder(elemn,av));
       sumy = sumy + xg(2,noder(elemn,av));
    nodenum = num2str(noder(elemn,av));
  text(xg(1,noder(elemn,av)),xg(2,noder(elemn,av)),nodenum,'Fontsize',9);
   end
   avercordx  = sumx/4;  % average x coordinate of quadrilateral
   avercordy  = sumy/4;  % average y coordinate of quadrilateral
  
   text(avercordx,avercordy,num2str(elemn),'Fontsize',7.5);




   

   
   end            
end
CONN = noder;
axis off

%% Force Vector
% Point Load and Moments
CORX = xg(1,:)';
CORY = xg(2,:)';
%% Assembly of Stiffness Matrix and Force Vector
SKE = NNE*NDOF; % Size of Element Stiffness Matrix
SKG = NN*NDOF; % Size of Global Stiffness Matrix
KG = zeros(SKG);
FG = zeros(SKG,1);
DOFM = zeros(NE,NDOF*NNE);
Ni = zeros(1,NNE);
dNizeta = zeros(1,NNE);
dNieta = zeros(1,NNE);
scf = 5/6;
for ne = 1:NE
FE = zeros(SKE,1);
Kemb = zeros(SKE,SKE);
Kegamma  = zeros(SKE,SKE);
  
E(ne) = Y;   
N1 = CONN(ne,1);
N2 = CONN(ne,2);
N3 = CONN(ne,3);
N4 = CONN(ne,4);
x1 = CORX(N1,1);
x2 = CORX(N2,1);
x3 = CORX(N3,1);
x4 = CORX(N4,1);
y1 = CORY(N1,1);
y2 = CORY(N2,1);
y3 = CORY(N3,1);
y4 = CORY(N4,1);

xn = [x1 x2 x3 x4];
yn = [y1 y2 y3 y4];
    

Dm = E(ne)*h/((1 - NU^2))*[1   NU   0 ;...
                                NU   1    0 ;...
                                0    0 (1-NU)/2];    
Db = E(ne)*h^3/(12*(1 - NU^2))*[1   NU   0 ;...
                                NU   1    0 ;...
                                0    0 (1-NU)/2];
Ds = scf*E(ne)*h/(1 - NU^2)*[(1-NU)/2 0; 0 (1-NU)/2];


for it = 1:2   
if it == 1
    ng = 3;
elseif it ==2
    ng = 1;
end               
[w, gp] = gausspoints(ng);
r = gp;
s = gp;


for ii = 1:ng
for jj = 1:ng
    
zeta = r(ii);
eta = s(jj);
weight = w(ii)*w(jj);    
    
    
zetai = [-1 1 1 -1];
etai  = [-1 -1 1 1];


for nc = 1:4   
    Ni(nc) = 1/4*(1 + zeta*zetai(nc))*(1 + eta*etai(nc));
    dNizeta(nc) = 1/4*zetai(nc)*(1 + eta*etai(nc));
    dNieta(nc)  = 1/4*etai(nc)*(1 + zeta*zetai(nc));       
end

%Jacobian Matrix
psixzeta = [dNizeta(1) 0 dNizeta(2) 0 dNizeta(3) 0 dNizeta(4) 0];
psiyzeta = [0 dNizeta(1) 0 dNizeta(2) 0 dNizeta(3) 0 dNizeta(4)];
psixeta = [dNieta(1) 0 dNieta(2) 0 dNieta(3) 0 dNieta(4) 0];
psiyeta = [0 dNieta(1) 0 dNieta(2) 0 dNieta(3) 0 dNieta(4)];
Xe = [xn(1) yn(1) xn(2) yn(2) xn(3) yn(3) xn(4) yn(4)];
J11 = psixzeta*Xe';
J12 = psiyzeta*Xe';
J21 = psixeta*Xe';
J22 = psiyeta*Xe';

J = [J11 J12;J21 J22];
jac = det(J);



 % Shape Functions     
psiu = [Ni(1) zeros(1,4) Ni(2) zeros(1,4) Ni(3) zeros(1,4) Ni(4) zeros(1,4)];
psiv  = [zeros(1,1) Ni(1) zeros(1,3) zeros(1,1) Ni(2) zeros(1,3) zeros(1,1) Ni(3) zeros(1,3) zeros(1,1) Ni(4) zeros(1,3)];
psiw   = [zeros(1,2) Ni(1) zeros(1,2) zeros(1,2) Ni(2) zeros(1,2) zeros(1,2) Ni(3) zeros(1,2) zeros(1,2) Ni(4) zeros(1,2)];
psithetax = [zeros(1,3) Ni(1) zeros(1,1) zeros(1,3) Ni(2) zeros(1,1) zeros(1,3) Ni(3) zeros(1,1) zeros(1,3) Ni(4) zeros(1,1)];
psithetay = [zeros(1,4) Ni(1) zeros(1,4) Ni(2) zeros(1,4) Ni(3) zeros(1,4) Ni(4)];
  
% Derivatives in zeta direction
dpsiuzeta = [dNizeta(1) zeros(1,4) dNizeta(2) zeros(1,4) dNizeta(3) zeros(1,4) dNizeta(4) zeros(1,4)];
dpsivzeta  = [zeros(1,1) dNizeta(1) zeros(1,3) zeros(1,1) dNizeta(2) zeros(1,3) zeros(1,1) dNizeta(3) zeros(1,3) zeros(1,1) dNizeta(4) zeros(1,3)];
dpsiwzeta  = [zeros(1,2) dNizeta(1) zeros(1,2) zeros(1,2) dNizeta(2) zeros(1,2) zeros(1,2) dNizeta(3) zeros(1,2) zeros(1,2) dNizeta(4) zeros(1,2)];
dpsithetaxzeta = [zeros(1,3) dNizeta(1) zeros(1,1) zeros(1,3) dNizeta(2) zeros(1,1) zeros(1,3) dNizeta(3) zeros(1,1) zeros(1,3) dNizeta(4) zeros(1,1)];
dpsithetayzeta= [zeros(1,4) dNizeta(1) zeros(1,4) dNizeta(2) zeros(1,4) dNizeta(3) zeros(1,4) dNizeta(4)];

% Derivatives in eta direction
dpsiueta = [dNieta(1) zeros(1,4) dNieta(2) zeros(1,4) dNieta(3) zeros(1,4) dNieta(4) zeros(1,4)];
dpsiveta  = [zeros(1,1) dNieta(1) zeros(1,3) zeros(1,1) dNieta(2) zeros(1,3) zeros(1,1) dNieta(3) zeros(1,3) zeros(1,1) dNieta(4) zeros(1,3)];
dpsiweta  = [zeros(1,2) dNieta(1) zeros(1,2) zeros(1,2) dNieta(2) zeros(1,2) zeros(1,2) dNieta(3) zeros(1,2) zeros(1,2) dNieta(4) zeros(1,2)];
dpsithetaxeta = [zeros(1,3) dNieta(1) zeros(1,1) zeros(1,3) dNieta(2) zeros(1,1) zeros(1,3) dNieta(3) zeros(1,1) zeros(1,3) dNieta(4) zeros(1,1)];
dpsithetayeta= [zeros(1,4) dNieta(1) zeros(1,4) dNieta(2) zeros(1,4) dNieta(3) zeros(1,4) dNieta(4)];




Bm = 1/jac*[J22 -J12 0 0; 0 0 -J21 J11; -J21 J11 J22 -J12]*[dpsiuzeta;dpsiueta;dpsivzeta ;dpsiveta];
Bb = 1/jac*[J22 -J12 0 0; 0 0 -J21 J11; -J21 J11 J22 -J12]*[dpsithetaxzeta;dpsithetaxeta;dpsithetayzeta ;dpsithetayeta];

Bgamma = inv(J)*[dpsiwzeta;dpsiweta] - [psithetax;psithetay];


 if it == 1
 fv = [zeros(2,1); Ni(1)*q0; zeros(2,1); zeros(2,1); Ni(2)*q0 ;zeros(2,1);zeros(2,1); Ni(3)*q0 ;zeros(2,1);zeros(2,1); Ni(4)*q0; zeros(2,1)]; 
 Kemb = Kemb + weight*jac*(Bm'*Dm*Bm +  Bb'*Db*Bb);
 FE = FE + jac*weight*fv;
 elseif it == 2
 Kegamma =  Kegamma + weight*jac*(Bgamma'*Ds*Bgamma);          
 end
end
end
end
KE = Kemb + Kegamma;
for j = 1:NNE
  nodenum = CONN(ne,j);
for k = 1:NDOF
 DOFM(ne,NDOF*(j-1) + k) = NDOF*(nodenum - 1) + k;
end
end

    [KG, FG] = planestressassemblyr(KE,FE,NNE,NDOF,DOFM(ne,:),KG,FG);
end
 
KGS = KG;
FGS =  FG;


%% Identify Nodes on Edges

% Bottom Edge
nodebottom = 1:gtrapy + 1;

%Top Edge
nodetop = nodebottom + (gtrapx + 1)*gtrapy; 

% Right Edge
noderight = gtrapx+1:gtrapx + 1:numnod;

% Left Edge
nodeleft = noderight - gtrapx; 

%% Application of Boundary Conditions
typebc =1;
%SS1
if typebc == 1
cn1 = [nodeleft noderight];
cdof1 = [5*cn1 - 3,5*cn1 - 2,5*cn1];
for ii = cdof1
KGS(ii,:) = 0;
KGS(:,ii) = 0;
KGS(ii,ii) = 1;
FGS(ii,1) = 0;
end

cn2 = [nodetop nodebottom];
cdof2 = [5*cn2 - 4,5*cn2 - 2,5*cn2-1];
for ii = cdof2
KGS(ii,:) = 0;
KGS(:,ii) = 0;
KGS(ii,ii) = 1;
FGS(ii,1) = 0;
end

elseif typebc == 3
%SS3
cn = [nodeleft noderight nodetop nodebottom];
cdof = [5*cn - 4,5*cn - 3,5*cn - 2];
for ii = cdof
KGS(ii,:) = 0;
KGS(:,ii) = 0;
KGS(ii,ii) = 1;
FGS(ii,1) = 0;
end
end

display ('Displacement Vector :')


u = KGS\FGS; % Displacement variables ;

dispw = u(3:5:end);

maxw = max(abs(dispw));

DD = Db(1,1);
nondw = maxw*DD/(q0*a^4)*1000
plotw = reshape(dispw,gtrapx+1,gtrapy+1);
figure(3); surf(plotw); xlabel("Plate x-dimension")
 ylabel("Plate y-dimension")
 zlabel("Displacement(mm)")

