function A=load_spiral_matrix(n)
% Loads a matrix with a spiral-like 
% spectrum. This is the discretization
% of an compact operator with a clustering
% point at l=1

tv=linspace(0,2*pi,n);
zv=3+((1./(1.5-tv/(2*pi))).^5).*(sin(5*tv)+0.5i*cos(5*tv))/30;
m=length(zv);
randn('seed',0);
A=randn(m);
A=A\diag(zv)*A;
