#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
        function TwoTensorField(field::Array)
          
        Constructor for the TwoTensorField struct. 

"""
     function TwoTensorField(field::Array{T}) where T <: Number

        el   = eltype(field)
        s    = collect(Int,size(field))
        @assert length(s) == 3
        
        dims = 2
        p    = s[1:2]
        nel  = s[3]

        return TwoTensorField{el}(dims,p,nel,field)
      end        

#---------------------------------------------------------------------- 
"""
        function NekField()
          
        Constructor for the Empty NekField struct. 

"""
      function NekField()

        T   = Float64
        hdr=""
        version=""
        wdsize=0
        lx1=0
        ly1=0
        lz1=0
        nel=0
        nelgt=0
        time=T(0)
        istep=0
        fid0=0
        nfileo=0
        rdcode=""
        p0th=T(0)
        ifprmesh=false
        glnum=Vector{Int}(undef,1)
        x=Array{T}(undef,0,0)
        y=Array{T}(undef,0,0)
        z=Array{T}(undef,0,0)
        u=Array{T}(undef,0,0)
        v=Array{T}(undef,0,0)
        w=Array{T}(undef,0,0)
        p=Array{T}(undef,0,0)
        t=Array{T}(undef,0,0)


       
        return NekField{T}(hdr,version,wdsize,lx1,ly1,lz1,nel,nelgt,time,
                          istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,
                          x,y,z,u,v,w,p,t)

     end        
#---------------------------------------------------------------------- 
"""
        function Re2Field()
          
        Constructor for the Empty Re2Field struct. 

"""
      function Re2Field()

        T   = Float64
        wdsize=0
        hdr=""
        version=""
        nelgt=0
        ldimr=0
        nelgv=0
        xc=Array{T}(undef,0,0)
        yc=Array{T}(undef,0,0)
        zc=Array{T}(undef,0,0)
        ncurve=0
        curveieg=Vector{Int}(undef,0)
        curveiside=Vector{Int}(undef,0)
        curveparam=Array{T}(undef,0,0)
        curvetype=Vector{String}(undef,0)
        cbl=Array{String}(undef,0,0)
        bl=Array{T}(undef,0,0)
       
        return Re2Field{T}(wdsize,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,
                           ncurve,curveieg,curveiside,curveparam,curvetype,
                           cbl,bl)
      end        
#---------------------------------------------------------------------- 
"""
        function ma2Hdr()
          
        Constructor for the Empty ma2Hdr struct. 

"""
      function ma2Hdr()
        
        return ma2Hdr("",0,0,0,0,0,0,0)
      end        
#---------------------------------------------------------------------- 
"""
        function ma2Field()
          
        Constructor for the Empty ma2Field struct. 

"""
      function ma2Field()
        
        return ma2Field(Vector{Int}(undef,0),Matrix{Int}(undef,0,0))
      end        
#---------------------------------------------------------------------- 




