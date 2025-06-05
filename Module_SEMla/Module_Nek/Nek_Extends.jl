# This is where we add extensions of Julia functions
#---------------------------------------------------------------------- 
function copy(f::NekFldHdr)

  T         = eltype(f.time)
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

  fnew = NekFldHdr{T}(version,wdsize,lx1,ly1,lz1,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh) 

  return fnew
end
#----------------------------------------------------------------------  

function copy(f::NekField)

  T         = eltype(f.hdr.time)
  hdr       = copy(f.hdr)
  glnum     = Base.copy(f.glnum)
  x         = Base.copy(f.x)
  y         = Base.copy(f.y)
  z         = Base.copy(f.z)
  u         = Base.copy(f.u)
  v         = Base.copy(f.v)
  w         = Base.copy(f.w)
  p         = Base.copy(f.p)
  t         = Base.copy(f.t)

  fnew = NekField{T}(hdr,glnum,x,y,z,u,v,w,p,t) 

  return fnew
end
#----------------------------------------------------------------------  
function copy(f::Re2Hdr)

  version   = f.version
  wdsize    = Base.copy(f.wdsize)
  nelgt     = Base.copy(f.nelgt)
  ldim      = Base.copy(f.ldim)
  nelgv     = Base.copy(f.nelgv)

  fnew      = Re2Hdr(version,wdsize,nelgt,ldim,nelgv)  

  return fnew
end
#----------------------------------------------------------------------
function copy(f::Re2Field)

  hdr             = copy(f.hdr)
  T               = eltype(f.xc)
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
  fnew            = Re2Field{T}(hdr,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl)  

  return fnew
end
#---------------------------------------------------------------------- 
function copy(f::Ma2Hdr)

  version   = f.version
  nel       = Base.copy(f.nel)
  nactive   = Base.copy(f.nactive)
  depth     = Base.copy(f.depth)
  d2        = Base.copy(f.d2)
  npts      = Base.copy(f.npts)
  nrank     = Base.copy(f.nrank)
  noutflow  = Base.copy(f.noutflow)

  fnew      = Ma2Hdr(version,nel,nactive,depth,d2,npts,nrank,noutflow)  

  return fnew
end
#----------------------------------------------------------------------
function copy(f::Ma2Field)

  hdr       = copy(f.hdr)
  pmap      = Base.copy(f.pmap)
  vmap      = Base.copy(f.vmap)

  fnew      = Ma2Field(hdr,pmap,vmap)  

  return fnew
end
#----------------------------------------------------------------------


