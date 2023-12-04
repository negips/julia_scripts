# This is where we add extensions of Julia functions

function copy(f::NekField)

  T         = eltype(f.time)
  hdr       = f.hdr
  version   = f.version
  wdsize    = Base.copy(f.wdsize)
  lx1       = Base.copy(f.lx1)
  ly1       = Base.copy(f.ly1)
  lz1       = Base.copy(f.lz1)
  nel       = Base.copy(f.nel)
  nelgt     = Base.copy(f.nelgt)
  time      = Base.copy(f.time)
  istep     = Base.copy(f.istep)
  fid0      = Base.copy(f.fid0)
  nfileo    = Base.copy(f.nfileo)
  rdcode    = Base.copy(f.rdcode)
  p0th      = Base.copy(f.p0th)
  ifprmesh  = Base.copy(f.ifprmesh)
  glnum     = Base.copy(f.glnum)
  x         = Base.copy(f.x)
  y         = Base.copy(f.y)
  z         = Base.copy(f.z)
  u         = Base.copy(f.u)
  v         = Base.copy(f.v)
  w         = Base.copy(f.w)
  p         = Base.copy(f.p)
  t         = Base.copy(f.t)

  fnew = NekField{T}(hdr,version,wdsize,lx1,ly1,lz1,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t) 

  return fnew
end
#----------------------------------------------------------------------  
function copy(f::Re2Field)

  T               = eltype(f.xc)
  wdsize          = Base.copy(f.wdsize)
  hdr             = f.hdr
  version         = f.version
  nelgt           = Base.copy(f.nelgt)
  ldimr           = Base.copy(f.ldimr)
  nelgv           = Base.copy(f.nelgv)
  xc              = Base.copy(f.xc)
  yc              = Base.copy(f.yc)
  zc              = Base.copy(f.zc)
  ncurve          = Base.copy(f.ncurve)            
  curveieg        = Base.copy(f.curveieg)
  curveiside      = Base.copy(f.curveiside)
  curveparam      = Base.copy(f.curveparam)
  curvetype       = Base.copy(f.curvetype)
  cbl             = Base.copy(f.cbl)
  bl              = Base.copy(f.bl)
  fnew            = Re2Field{T}(wdsize,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl)  

  return fnew
end
#---------------------------------------------------------------------- 
function copy(f::ma2Field)

  pmap      = Base.copy(f.pmap)
  vmap      = Base.copy(f.vmap)

  fnew      = ma2Field(pmap,vmap)  

  return fnew
end
#----------------------------------------------------------------------

# Re2Field{Float32}(wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl)

