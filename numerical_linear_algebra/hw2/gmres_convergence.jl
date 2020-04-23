#!/usr/bin/julia

function gmres_convergence(A)

  evals = eigvals(Matrix(A));
  
  er = real(evals);
  ei = imag(evals);
  
  rmax,imax = findmax(er);
  rmin,imin = findmin(er);
  
# Creating a more sophisticated circle
  ind = sortperm(er);
  ers = er[ind];
  eis = ei[ind];
  neig = length(ers);
  rmax = ers[neig-1];
  rmin = ers[1];
  l1 = rmax - rmin;
  l2 = abs(eis[2]);
  tanb = l2/l1;
  l3 = l1*tanb^2;
  dia = l1 + l3;
  rad = dia/2;
  opp_x = rmax - dia;
  cen = (rmax + opp_x)/2;
  cfac = rad/cen;

  c2 = ers[neig];
  
  rho = sqrt(rad*(rad + abs(cen-c2))/abs(cen*c2))

#  conv2 = 10000. *r[1]*(cfac).^ik;
#  conv3 = 1000. *r[1]*(rho).^(ik.-1);

  return cfac,rad,cen

end






