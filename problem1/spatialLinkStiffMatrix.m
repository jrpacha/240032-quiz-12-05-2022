function Ke=spatialLinkStiffMatrix(nod,elem,e,E,A)
%
% Matriu de rigidesa per un element Barra en 2-dim
% 
% (c) Numerical factory
%
x1=nod(elem(e,1),1);
y1=nod(elem(e,1),2);
z1=nod(elem(e,1),3);
x2=nod(elem(e,2),1);
y2=nod(elem(e,2),2);
z2=nod(elem(e,2),3);

x21=x2-x1;
y21=y2-y1;
z21=z2-z1;

Le = sqrt(x21*x21+y21*y21+z21*z21);
coef = ((E(e)*A(e))/Le^3);
M = [x21*x21, x21*y21, x21*z21;
y21*x21, y21*y21, y21*z21;
z21*x21, z21*y21, z21*z21];
Ke = [ M, -M;
      -M, M];
Ke = coef*Ke;
