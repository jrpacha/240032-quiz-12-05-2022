clearvars
close all

% ====================== Quiz 12-05-2022 Problema 2 =======================
% ================================ Data ===================================
kc = 1.0;            %Thermal conductivity
f = 1.0;
tempLeftVert=133.0;  %ºC. Temperature on the left vertical boundary (LVB)
tempRightHoriz=50.0; %ºC. Temperature on the right horizontal boudary (RHB) 

q0 = -0.40;          %For part (c): Constant flow per unit length on 
                     %the circle exterior boundary (CEB)

% ================================ Mesh ===================================
fileSols = 'solutionP2.txt'; 

eval('meshClipQuiz');
[nodes,elem] = generateTriangFromQuadMesh(nodes,elem);
numNodes = size(nodes,1);
numElem = size(elem,1);

% Plot the mesh
numbering= 0; %= 1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering);


%%% Parts (a) and (b)

% Select Boundary points
indNodLV = find(nodes(:,1) < 0.01); %indices of nodes at LVB
indNodRH = find(nodes(:,2) < 0.01); %indices of nodes at RHB

hold on %plot the points corresponding to the selected indexes
title('Selected boundary nodes for parts (a) and (b)')
plot(nodes(indNodLV,1),nodes(indNodLV,2),'ok','lineWidth',1,...
    'markerFaceColor','red','markerSize',5)
plot(nodes(indNodRH,1),nodes(indNodRH,2),'ok','lineWidth',1,...
    'markerFaceColor','blue','markerSize',5)
hold off

%Define the coefficients vector of the model equation
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNodes);    %global stiff matrix
F=zeros(numNodes,1);  %global internal forces vector
Q=zeros(numNodes,1);  %global secondary variables vector
u=zeros(numNodes,1);  %solutions

for e = 1:numElem
    [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e);
    rows= [elem(e,1); elem(e,2); elem(e,3)];
    cols= rows;
    K(rows,cols)= K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)= F(rows) + Fe;
    end
end

%Booundary Conditions
fixedNods = [indNodLV', indNodRH']; %fixed Nodes (global numbering)
%fixedNods = unique(fixedNods);
freeNods = setdiff(1:numNodes,fixedNods); %free Nodes (global numbering)

% Essential B.C.
u(indNodLV) = tempLeftVert;   %temperature on the LVB
u(indNodRH) = tempRightHoriz; %temperature on the RVB

%Reduced system
Fm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

%Compute the solution
um = Km\Fm;
u(freeNods)= um;

%PostProcess: Compute secondary variables and plot the temperature
%distribution
Q = K*u - F;

titol='Temperature Distribution parts (a), (b)';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

%%% Solutions of parts (a) and (b)
clc
%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
fOut=fopen(fileSols,'w');
fprintf(fOut,'\t\t   *** Problem 2 ***\n');
fprintf(fOut,...
    '(a) The entry (400,400) of the global stiffness matrix is\n');
fprintf(fOut, ...
    '    K(400,400) = %.4e\n\n',K(400,400));
fprintf(fOut, ...
    '    Hint1. The value of K(300,300) = %.6e\n\n',K(300,300));
fprintf(fOut, ...
    '(b) The minimum of the values of Q on the left vertical boundary\n');
fprintf(fOut, ...
    '    is %.4e\n\n',min(Q(indNodLV)));
fprintf(fOut, ...
    '    Hint2. The maximum in absolute value of Q over all the nodes\n');
fprintf(fOut, ...
    '    is %.4e\n\n',max(abs(Q)));

%format short e
%K(400,400)
%K(300,300)
%min(Q(indNodLV))
%max(abs(Q))

%%% Part (c)
indNodCirc = find(sqrt(nodes(:,1).^2+nodes(:,2).^2) > 46.01); %Nods at CEB

figure() %plot the points corresponding to the selected indexes
plotElementsOld(nodes,elem,numbering);
title('Selected boundary nodes for part (c)')
hold on 
plot(nodes(indNodLV,1),nodes(indNodLV,2),'ok','lineWidth',1,...
    'markerFaceColor','red','markerSize',5)
plot(nodes(indNodRH,1),nodes(indNodRH,2),'ok','lineWidth',1,...
    'markerFaceColor','blue','markerSize',5)
plot(nodes(indNodCirc,1),nodes(indNodCirc,2),'ok','lineWidth',1,...
    'markerFaceColor','green','markerSize',5)
hold off

clear Q;
Q = zeros(numNodes,1); %Reset vector Q

%Natural B.C.
indBC = indNodCirc';
Q=applyConstantNaturalBC(nodes,elem,indBC,q0,Q);

%Essential B.C.
u(indNodLV) = tempLeftVert;
u(indNodRH) = tempRightHoriz;

%Reduced system
Fm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

%Compute the solution
um = Km\Fm;
u(freeNods)= um;

%PostProcess: Compute secondary variables and plot the temperature
%distribution
Q = K*u - F;

titol='Temperature Distribution (c)';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

%%% Solution of part (c)
%--------------------------------------------------------------------------
%     Fancy output. Do not waste your time with this at the exams!
%--------------------------------------------------------------------------
fprintf(fOut, ...
    '(c) Now consider that the circle exterior boundary is no more\n');
fprintf(fOut, ...
    '   isolated, but there is a constant flow per unit length qn = \n');
fprintf(fOut, ...
    '   %.2f trought that boundary\n\n',q0);
fprintf(fOut, ...
    '   Then, the minimum of the values of Q on the left vertical\n');
fprintf(fOut, ...
    '   boundary is: %.4e\n\n',min(Q(indNodLV)));
fprintf(fOut, ...
    '   Hint3. The maximum in absolute value of Q over all the nodes\n');
fprintf(fOut, ...
    '   is %.4e\n',max(abs(Q)));
% -------------------------------------------------------------------------
%min(Q(indNodLV))
%max(abs(Q))

fclose(fOut);
type(fileSols);

