%%
clearvars
close all

% ===================== Quiz 12-05-2022 Problem 1 =========================
fileSols = 'solutionP1.txt';

%Real constants: Materials and sections area
Area = 250.0;            %section Area (in mm^2);
newArea = 150.0;         %in part (c), the section area of all the bars
                         %except the ones at the basis is changed to 
                         % 150.0 mm^2
Y=1.5e5;                 %Young modulus of the bars, in N/mm^2

%Goemetry
ll = 2755;               %mm
h = 3414;                %mm

H=ll*sqrt(3.0)/2.0;      %mm
D=ll/2;                  %mm

%F=[29705;-20515;0.0]; %force applied to node 6
F=[30000; 30000; -30000];

%Output Data
fOut=fopen(fileSols,'w');
clc
%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
fprintf(fOut,'\t\t    ** Problem 1 **\n');
fprintf(fOut,'Length of the bars at the basis and at the top:\n');
fprintf(fOut,'       l = %e mm\n',ll);
fprintf(fOut,'Length of the vertical bars:\n');
fprintf(fOut,'       h = %e mm\n',h);
fprintf(fOut,'Section area of the bars:\n');
fprintf(fOut,'       A = %e mm^2\n',Area);
fprintf(fOut,'Young modulus:\n');
fprintf(fOut,'       Y = %e N/mm^2\n',Y);
fprintf(fOut,'Force at node 5:\n');
fprintf(fOut,'       F = [%e, %e, %e]\n',F);
fprintf(fOut,'\n');
%--------------------------------------------------------------------------

% nodes=[0.0,0.0,0.0;
%        3600.0,0.0,0.0;
%        7200.0,0.0,0.0;
%        10800.0,0.0,0.0;
%        1800.0,H,0.0;
%        5400.0,H,z;
%        9000.0,H,0.0];

nodes = [
    0, 0, 0;
    ll, 0, 0;
    D, H, 0;
    0, 0, h;
    ll, 0, h;
    D, H, h];

elem = [
    1, 2;
    2, 3;
    1, 3;
    1, 4;
    2, 5;
    3, 6;
    4, 5;
    5, 6;
    4, 6];
      
numNodes=size(nodes,1);
numElem=size(elem,1);
ndim=size(nodes,2);

numbering = 1;
%plotElements(nodes, elem, numbering);
plotElementsOld(nodes, elem, numbering);

%Real constants
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(ndim*numNodes,1);
Q=zeros(ndim*numNodes,1);
K=zeros(ndim*numNodes);

%%% Part (A)
%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
Ke=spatialLinkStiffMatrix(nodes,elem,3,E,A);
fprintf(fOut,...
    '(a) The value of the component (4,2) of the local stiffness\n');
fprintf(fOut,...
    '    matrix for the 3rd element, K3, is K3(4,2) = %.4e\n\n',Ke(4,2));
fprintf(fOut,...
    '    Hint1. The value of K3(2,2) = %.6e\n\n',Ke(2,2));
%--------------------------------------------------------------------------

%%% Part (B)
for e=1:numElem
    Ke=spatialLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-2,ndim*elem(e,1)-1,ndim*elem(e,1),...
          ndim*elem(e,2)-2,ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Boundary Conditions
nod = 5;   %Only one free node
freeNods = [ndim*nod-2,ndim*nod-1,ndim*nod]; 
fixedNods = setdiff(1:ndim*numNodes,freeNods);

%Natural B.C. (Loads)
Q(ndim*nod-2)=F(1); %N
Q(ndim*nod-1)=F(2); %N
Q(ndim*nod)=F(3);   %N 

%Essential B.C
u(fixedNods) = 0; %redundant, since the vector has been already
                  %initialised to 0

%Reduced system: Km * um = Qm
Qm=Q(freeNods)-K(freeNods,fixedNods)*u(fixedNods); 
   %Remark: since in this case u(fixedNods) = 0, the term
   %        -K(freeNods,fixedNods)*u(fixedNods) is no necessary and and 
   %        we could have written Qm = Q(freeNods)   
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

U = [u(1:ndim:end),u(2:ndim:end),u(3:ndim:end)];

%Post porcess: plot the deformed structure
esc=20;
plotDeformedTruss(nodes, elem, u, esc);
title('Deformed structure parts (a) and (b)')
hold off

%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
fprintf(fOut,...
    '(b) The absolute value of the y-displacement of the node 5\n');
fprintf(fOut,...
    '    is: |UY(5)| = %.4e\n\n',abs(U(5,2)));
fprintf(fOut,...
    '    Hint2. The absolute value of the z-displacement of the\n');
fprintf(fOut,...
    '    node 5 is |UZ(5)| = %.4e\n\n',abs(U(5,3)));
%--------------------------------------------------------------------------

%%% Part (C)
% Execp the three bars at the basis, the section area of all the other 
% is changed to 150mm^2

A = newArea*ones(1,numElem);
A(1:3) = Area;

clear K;
K = zeros(ndim*numNodes); %Reset matrix K to zero

for e=1:numElem
    Ke=spatialLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-2,ndim*elem(e,1)-1,ndim*elem(e,1),...
          ndim*elem(e,2)-2,ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

Km = K(freeNods,freeNods); %Remark: neither u(freeNods) nor Qm change!
um=Km\Qm;
u(freeNods) = um;

U = [u(1:ndim:end),u(2:ndim:end),u(3:ndim:end)];

%Post porcess: plot the deformed structure
esc=20;
plotDeformedTruss(nodes, elem, u, esc);
title('Deformed structure part (c)')
hold off

%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
fprintf(fOut,...
    '(c) Now, except the three bars at the bassis, the section\n');
fprintf(fOut,...
    '    area of all the others is changed to A = %.1f mm^2.\n',newArea);
fprintf(fOut,'    Then:\n\n');
fprintf(fOut,...
    '    The absolute value of the y-displacement of the node 5\n');
fprintf(fOut,...
    '    is |UY(5)| = %.4e\n\n',abs(U(5,2)));
fprintf(fOut,...
    '    Hint3. The absolute value of the x-displacement of the\n');
fprintf(fOut,...
    '    node 5 is |UX(5)| = %.4e\n',abs(U(5,1)));
%--------------------------------------------------------------------------
fclose(fOut);
type(fileSols);