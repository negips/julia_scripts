#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
        function NekFldHdr()
          
        Constructor for the Empty NekFldHdr struct. 

"""
      function NekFldHdr(T::DataType)

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
       
        return NekFldHdr{T}(version,wdsize,lx1,ly1,lz1,nel,nelgt,time,
                          istep,fid0,nfileo,rdcode,p0th,ifprmesh)

     end        
#---------------------------------------------------------------------- 

"""
        function NekField()
          
        Constructor for the Empty NekField struct. 

"""
      function NekField(T::DataType)

        hdr       = NekFldHdr(T)
        glnum     = Vector{Int}(undef,0)
        x         = Array{T}(undef,0,0)
        y         = Array{T}(undef,0,0)
        z         = Array{T}(undef,0,0)
        u         = Array{T}(undef,0,0)
        v         = Array{T}(undef,0,0)
        w         = Array{T}(undef,0,0)
        p         = Array{T}(undef,0,0)
        t         = Array{T}(undef,0,0)
       
        return NekField{T}(hdr,glnum,x,y,z,u,v,w,p,t)

     end        
#---------------------------------------------------------------------- 
"""
        function Re2Hdr()
          
        Constructor for the Empty Re2Hdr struct. 

"""
      function Re2Hdr()

        version=""
        wdsize=0
        nelgt=0
        ldim=0
        nelgv=0
       
        return Re2Hdr(version,wdsize,nelgt,ldim,nelgv)
      end        
#---------------------------------------------------------------------- 

"""
        function Re2Field()
          
        Constructor for the Empty Re2Field struct. 

"""
      function Re2Field()

        T   = Float64
        hdr = Re2Hdr()
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
       
        return Re2Field{T}(hdr,xc,yc,zc,ncurve,curveieg,curveiside,
                           curveparam,curvetype,cbl,bl)

      end        
#---------------------------------------------------------------------- 

"""
        function Ma2Hdr()
          
        Constructor for the Empty ma2Hdr struct. 

"""
      function Ma2Hdr()
        
        return Ma2Hdr("",0,0,0,0,0,0,0)
      end        
#---------------------------------------------------------------------- 
"""
        function Ma2Field()
          
        Constructor for the Empty ma2Field struct. 

"""
      function Ma2Field()
        
        return Ma2Field(Ma2Hdr(),Vector{Int}(undef,0),Matrix{Int}(undef,0,0))
      end        
#---------------------------------------------------------------------- 

"""
        function ReaData()
          
        Constructor for the Empty ReaData struct. 

"""
      function ReaData()

        T               = Float64
        re2             = Re2Field()
        npars           = 0
        par_vals        = zeros(Float64,npars)
        par_names       = Vector{String}(undef,npars)
        par_desc        = Vector{String}(undef,npars)
        nrst            = 0
        rstfiles        = Vector{String}(undef,nrst)
        rstoptions      = Vector{String}(undef,nrst)      
        nic             = 0
        icfiles         = Vector{String}(undef,nic)
        ndf             = 0
        driveforce      = Vector{String}(undef,ndf)
        nhist           = 0
        history         = Vector{String}(undef,nhist)
        noutspec        = 0
        OutFlds         = Vector{String}(undef,noutspec)
        OutSpec         = Vector{Integer}(undef,noutspec)
        Objs            = NekObjects()
       
        return ReaData{T}(re2,par_vals,par_names,par_desc,nrst,rstfiles,rstoptions,nic,icfiles,ndf,driveforce,nhist,history,noutspec,OutFlds,OutSpec,Objs)

      end        
#---------------------------------------------------------------------- 
"""
        function NekObjects()
          
        Constructor for the Empty NekObjects struct. 

"""
      function NekObjects()
       
        return NekObjects(0,0,0,0)

      end        
#---------------------------------------------------------------------- 

#----------------------------------------------------------------------
# Non-Empty Constructors
#---------------------------------------------------------------------- 
"""
        function Re2Field(xc::Matrix{T},yc::Matrix{T},zc::Matrix{T},BC::Matrix{String},Par::Array{T}) where {T <: AbstractFloat}
          
        Constructor for Re2Field struct from Coordinates and Boundary Data 

"""
      function Re2Field(xc::Matrix{T},yc::Matrix{T},zc::Matrix{T},BC::Matrix{String},Par::Array{T}) where {T <: AbstractFloat}

        if T == Float64
          wdsize = 8
        else
          wdsize = 4
        end

        nc,nelgt = size(xc)
        if nc == 4
          ndim = 2
        elseif nc == 8
          ndim = 3
        else
          println("Unknown dim")
        end

        nelgv = nelgt
        version = "#v002"     # Standard

        re2hdr = Re2Hdr(version,wdsize,nelgt,ndim,nelgv)
        ncurve=0
        curveieg=Vector{Int}(undef,0)
        curveiside=Vector{Int}(undef,0)
        curveparam=Array{T}(undef,0,0)
        curvetype=Vector{String}(undef,0)
       
        return Re2Field{T}(re2hdr,xc,yc,zc,ncurve,curveieg,curveiside,
                           curveparam,curvetype,BC,Par)
      end        
#---------------------------------------------------------------------- 

"""
        function ReaData(re2::ReField,npars::Int64)
          
        Const ReaData struct from re2 data. 

"""
      function ReaData(re2::Re2Field{T},npars::Int64) where {T}

        pnam,pdesc      = Nek.param_descriptions()
        pn              = copy(pnam[1:npars])
        pd              = copy(pdesc[1:npars])
        pv              = zeros(Float64,npars)
       
        nrst            = 0
        rstfiles        = Vector{String}(undef,nrst)
        rstoptions      = Vector{String}(undef,nrst)      
        nic             = 0
        icfiles         = Vector{String}(undef,nic)
        ndf             = 0
        driveforce      = Vector{String}(undef,ndf)
        nhist           = 0
        history         = Vector{String}(undef,nhist)
        noutspec        = 0
        OutFlds         = Vector{String}(undef,noutspec)
        OutSpec         = Vector{Integer}(undef,noutspec)
        Objs            = NekObjects()
       
        return ReaData{T}(re2,pv,pn,pd,nrst,rstfiles,rstoptions,nic,icfiles,ndf,driveforce,nhist,history,noutspec,OutFlds,OutSpec,Objs)

      end        
#---------------------------------------------------------------------- 




