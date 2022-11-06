function ConsistentPoisson2D(vB,pB,binv,gradnodes,gradweights,divnodes,divweights,prec)

## Consistent Integration (2D)
#--------------------------------------------------

lxv = PolynomialBases.degree(vB) + 1
lxp = PolynomialBases.degree(pB) + 1
ng  = length(gradnodes)
nd  = length(divnodes)

dxV_g = zeros(prec,ng,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_g[:,j] = derivative_at(gradnodes,v,vB.nodes,vB.baryweights)
end
dxV_gT = transpose(dxV_g)

intV_g = zeros(prec,ng,lxv)
intV_g = interpolation_matrix(gradnodes,vB.nodes,vB.baryweights)
intV_g_2d  = kron(intV_g,intV_g);

intP_g = zeros(prec,ng,lxp)
intP_g = interpolation_matrix(gradnodes,pB.nodes,pB.baryweights)
intP_g_2d  = kron(intP_g,intP_g);

# Gradients
gw2d = kron(gradweights,gradweights)
wg2d_x = kron(intV_g',dxV_gT)*Diagonal(gw2d[:])*intP_g_2d
wg2d_y = kron(dxV_gT,intV_g')*Diagonal(gw2d[:])*intP_g_2d


dxV_d = zeros(prec,nd,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_d[:,j] = derivative_at(divnodes,v,vB.nodes,vB.baryweights)
end
dxV_dT = transpose(dxV_d)

intV_d = zeros(prec,nd,lxv)
intV_d = interpolation_matrix(divnodes,vB.nodes,vB.baryweights)
intV_d_2d  = kron(intV_d,intV_d);


intP_d = zeros(prec,nd,lxp)
intP_d = interpolation_matrix(divnodes,pB.nodes,pB.baryweights)
intP_d_2d  = kron(intP_d,intP_d);


dw2d = kron(divweights,divweights)
div2d_x  = (intP_d_2d')*Diagonal(dw2d)*kron(intV_d,dxV_d)*Diagonal(binv[:])
div2d_y  = (intP_d_2d')*Diagonal(dw2d)*kron(dxV_d,intV_d)*Diagonal(binv[:])

CPoisson2D = div2d_x*wg2d_x + div2d_y*wg2d_y
mass       = kron(pB.weights,pB.weights)
Mass       = Diagonal(mass[:]);

return CPoisson2D,Mass
end

#---------------------------------------------------------------------- 

function ConsistentPoisson2D_ρ(vB,pB,binv,ρ,gradnodes,gradweights,divnodes,divweights,prec)

## Consistent Integration (2D)
#--------------------------------------------------

lxv = PolynomialBases.degree(vB) + 1
lxp = PolynomialBases.degree(pB) + 1
ng  = length(gradnodes)
nd  = length(divnodes)

dxV_g = zeros(prec,ng,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_g[:,j] = derivative_at(gradnodes,v,vB.nodes,vB.baryweights)
end
dxV_gT = transpose(dxV_g)

intV_g = zeros(prec,ng,lxv)
intV_g = interpolation_matrix(gradnodes,vB.nodes,vB.baryweights)
intV_g_2d  = kron(intV_g,intV_g);

intP_g = zeros(prec,ng,lxp)
intP_g = interpolation_matrix(gradnodes,pB.nodes,pB.baryweights)
intP_g_2d  = kron(intP_g,intP_g);

# Gradients
gw2d = kron(gradweights,gradweights)
wg2d_x = kron(intV_g',dxV_gT)*Diagonal(gw2d[:])*intP_g_2d
wg2d_y = kron(dxV_gT,intV_g')*Diagonal(gw2d[:])*intP_g_2d


dxV_d = zeros(prec,nd,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_d[:,j] = derivative_at(divnodes,v,vB.nodes,vB.baryweights)
end
dxV_dT = transpose(dxV_d)

intV_d = zeros(prec,nd,lxv)
intV_d = interpolation_matrix(divnodes,vB.nodes,vB.baryweights)
intV_d_2d  = kron(intV_d,intV_d);


intP_d = zeros(prec,nd,lxp)
intP_d = interpolation_matrix(divnodes,pB.nodes,pB.baryweights)
intP_d_2d  = kron(intP_d,intP_d);

intVP_d     = zeros(prec,nd,lxp)
intVP_d     = interpolation_matrix(divnodes,vB.nodes,vB.baryweights)
intVP_d_2d  = kron(intVP_d,intVP_d);

ρ_d_2d      = intVP_d_2d*ρ[:]

dw2d = kron(divweights,divweights)
div2d_x  = (intP_d_2d')*Diagonal(ρ_d_2d)*Diagonal(dw2d)*kron(intV_d,dxV_d)*Diagonal(binv[:])*Diagonal(1.0./ρ[:])
div2d_y  = (intP_d_2d')*Diagonal(ρ_d_2d)*Diagonal(dw2d)*kron(dxV_d,intV_d)*Diagonal(binv[:])*Diagonal(1.0./ρ[:])

CPoisson2D = (div2d_x*wg2d_x + div2d_y*wg2d_y)
mass       = kron(pB.weights,pB.weights)
Mass       = Diagonal(mass[:]);

return CPoisson2D,Mass
end

#---------------------------------------------------------------------- 

function ConsistentPoisson2D_Dealias(vB,pB,binv,gradnodes,gradweights,divnodes,divweights,dealnodes,dealweights,prec)

## Consistent Integration (2D)
#--------------------------------------------------

lxv = PolynomialBases.degree(vB) + 1
lxp = PolynomialBases.degree(pB) + 1
ng  = length(gradnodes)
nd  = length(divnodes)
ndeal = length(dealnodes)     # Dealiasing nodes


# Gradients
dxV_g = zeros(prec,ng,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_g[:,j] = derivative_at(gradnodes,v,vB.nodes,vB.baryweights)
end
dxV_gT = transpose(dxV_g)

intV_g = zeros(prec,ng,lxv)
intV_g = interpolation_matrix(gradnodes,vB.nodes,vB.baryweights)
intV_g_2d  = kron(intV_g,intV_g);

intP_g = zeros(prec,ng,lxp)
intP_g = interpolation_matrix(gradnodes,pB.nodes,pB.baryweights)
intP_g_2d  = kron(intP_g,intP_g);

gw2d = kron(gradweights,gradweights)
wg2d_x = kron(intV_g',dxV_gT)*Diagonal(gw2d[:])*intP_g_2d
wg2d_y = kron(dxV_gT,intV_g')*Diagonal(gw2d[:])*intP_g_2d


# Dealiasing
m2d       = kron(vB.weights,vB.weights)
intV_deal = zeros(prec,ndeal,lxv)
intV_deal = interpolation_matrix(dealnodes,vB.nodes,vB.baryweights)
intV_deal_2d  = kron(intV_deal,intV_deal);

m2d_deal  = kron(dealweights,dealweights)
binv_deal = m2d_deal[:].*(intV_deal_2d*binv[:])
deal      = Diagonal(1.0./m2d[:])*(intV_deal_2d')*Diagonal(binv_deal)*intV_deal_2d


# Divergence
dxV_d = zeros(prec,nd,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_d[:,j] = derivative_at(divnodes,v,vB.nodes,vB.baryweights)
end
dxV_dT = transpose(dxV_d)

intV_d = zeros(prec,nd,lxv)
intV_d = interpolation_matrix(divnodes,vB.nodes,vB.baryweights)
intV_d_2d  = kron(intV_d,intV_d);


intP_d = zeros(prec,nd,lxp)
intP_d = interpolation_matrix(divnodes,pB.nodes,pB.baryweights)
intP_d_2d  = kron(intP_d,intP_d);


dw2d = kron(divweights,divweights)
div2d_x  = (intP_d_2d')*Diagonal(dw2d)*kron(intV_d,dxV_d)*deal
div2d_y  = (intP_d_2d')*Diagonal(dw2d)*kron(dxV_d,intV_d)*deal

CPoisson2D = div2d_x*wg2d_x + div2d_y*wg2d_y
mass       = kron(pB.weights,pB.weights)
Mass       = Diagonal(mass[:]);

return CPoisson2D,Mass
end

#---------------------------------------------------------------------- 

function ConsistentPoisson2D_Dealias_ρ(vB,pB,binv,ρ,gradnodes,gradweights,divnodes,divweights,dealnodes,dealweights,prec)

## Consistent Integration (2D)
#--------------------------------------------------

lxv = PolynomialBases.degree(vB) + 1
lxp = PolynomialBases.degree(pB) + 1
ng  = length(gradnodes)
nd  = length(divnodes)
ndeal = length(dealnodes)     # Dealiasing nodes


# Gradients
dxV_g = zeros(prec,ng,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_g[:,j] = derivative_at(gradnodes,v,vB.nodes,vB.baryweights)
end
dxV_gT = transpose(dxV_g)

intV_g = zeros(prec,ng,lxv)
intV_g = interpolation_matrix(gradnodes,vB.nodes,vB.baryweights)
intV_g_2d  = kron(intV_g,intV_g);

intP_g = zeros(prec,ng,lxp)
intP_g = interpolation_matrix(gradnodes,pB.nodes,pB.baryweights)
intP_g_2d  = kron(intP_g,intP_g);

gw2d = kron(gradweights,gradweights)
wg2d_x = kron(intV_g',dxV_gT)*Diagonal(gw2d[:])*intP_g_2d
wg2d_y = kron(dxV_gT,intV_g')*Diagonal(gw2d[:])*intP_g_2d


# Dealiasing
m2d       = kron(vB.weights,vB.weights)
intV_deal = zeros(prec,ndeal,lxv)
intV_deal = interpolation_matrix(dealnodes,vB.nodes,vB.baryweights)
intV_deal_2d  = kron(intV_deal,intV_deal);

m2d_deal    = kron(dealweights,dealweights)
binv_fine   = (intV_deal_2d*binv[:])
rhoinv_fine = (intV_deal_2d*(1.0./ρ[:]))
bρinv_deal  = m2d_deal[:].*(binv_fine.*rhoinv_fine)

deal      = Diagonal(1.0./m2d[:])*(intV_deal_2d')*Diagonal(bρinv_deal)*intV_deal_2d


# Divergence
dxV_d = zeros(prec,nd,lxv)
for j in 1:lxv
  v = zeros(prec,lxv)
  v[j] = 1.0
  dxV_d[:,j] = derivative_at(divnodes,v,vB.nodes,vB.baryweights)
end
dxV_dT = transpose(dxV_d)

intV_d = zeros(prec,nd,lxv)
intV_d = interpolation_matrix(divnodes,vB.nodes,vB.baryweights)
intV_d_2d  = kron(intV_d,intV_d);


intP_d = zeros(prec,nd,lxp)
intP_d = interpolation_matrix(divnodes,pB.nodes,pB.baryweights)
intP_d_2d  = kron(intP_d,intP_d);


dw2d = kron(divweights,divweights)
div2d_x  = (intP_d_2d')*Diagonal(dw2d)*kron(intV_d,dxV_d)*deal
div2d_y  = (intP_d_2d')*Diagonal(dw2d)*kron(dxV_d,intV_d)*deal

CPoisson2D = div2d_x*wg2d_x + div2d_y*wg2d_y
mass       = kron(pB.weights,pB.weights)
Mass       = Diagonal(mass[:]);

return CPoisson2D,Mass
end








