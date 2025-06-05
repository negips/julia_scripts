#     Author:     Prabal Negi
#     Nek based routines in SEMla
#---------------------------------------------------------------------- 
function write_reafile(case::String,rea::ReaData)

      # Open File
      fname = case*".rea"
      io    = open(fname, "w")

      ndim  = rea.re2.hdr.ldim
      ver   = rea.re2.hdr.version

      # Write Header
      nekver = 19.0
      write_nekhdr(io,nekver,ndim)

      # Parameter Values
      vals    = rea.par_vals
      names   = rea.par_names
      desc    = rea.par_desc
      write_params(io,vals,names,desc)

      # Passive Scalars
      nps = 0
      write_ps_hdr(io,nps)

      # Logical Switches
      nlogic = 0
      write_logicalswitches_hdr(io,nlogic)

      # Prenek Header
      write_prenek_hdr(io)

      # Mesh data
      write_mesh(io,rea.re2)

      # Restarts
      nrst        = rea.nrst
      rstfiles    = rea.rstfiles
      rstoptions  = rea.rstoptions
      write_restart(io,nrst,rstfiles,rstoptions)

      # Initial Condition
      nic         = rea.nic
      icfiles     = rea.icfiles
      write_ic(io,nic,icfiles)

      # Driving Force
      ndf         = rea.ndf
      driveforce  = rea.driveforce
      write_driveforce(io,ndf,driveforce)


      # No Variable property data
      nvp         = 1
      npackets    = 0
      datapacket  = Vector{String}(undef,npackets)
      write_varprop(io,nvp,npackets,datapacket)

      # History
      nhist       = rea.nhist
      history     = rea.history
      write_hist(io,nhist,history)

      # I/O Specification
      nio         = rea.noutspec
      outflds     = rea.outflds
      outspec     = rea.outspec
      write_outputspec(io,nio,outflds,outspec)

      # Objects
      objects     = rea.objects
      write_objectspec(io,objects)

      # Close file
      close(io)

      return nothing
end  
#---------------------------------------------------------------------- 
function write_mesh(fid::IOStream,re2::Re2Field)

      # Mesh
      # Header
      ifre2 = false
      ndim  = re2.hdr.ldim
      nelg  = re2.hdr.nelgt
      nelv  = re2.hdr.nelgv
      write_mesh_hdr(fid,ndim,nelg,nelv,ifre2)

      # Coordinates
      xc    = re2.xc
      yc    = re2.yc
      zc    = re2.zc
      write_mesh_coords(fid,ndim,xc,yc,zc)

      # Curved Side data
      ncurve      = re2.ncurve
      curveieg    = re2.curveieg
      curveedge   = re2.curveiside
      curveparams = re2.curveparam
      curvetype   = re2.curvetype

      # Header
      write_curve_hdr(fid,ncurve)

      # Curve Params
      write_curvedata(fid,nelg,ncurve,curveieg,curveedge,curveparams,curvetype)

      # Fluid BCs
      # Header
      write_fluidBC_hdr(fid)

      # BC data
      BC          = re2.cbl
      Pars        = re2.bl
      write_fluidBC(fid,ndim,nelg,BC,Pars)

      # Thermal BC
      write_zerothermalBC_hdr(fid)

      return nothing
end  
#---------------------------------------------------------------------- 
function write_nekhdr(fid::IOStream,ver::Float64,ndim::Int)

      hdr="****** PARAMETERS ******"
      @printf(fid,"%s\n",hdr);
      
      space5=blanks(5);
      hdr="NEKTON VERSION";
      @printf(fid,"%12.6f%s%s\n",ver,space5,hdr);

      space2=blanks(2);
      hdr=" DIMENSIONAL RUN";
      @printf(fid,"%s%2i%s\n",space2,ndim,hdr);

      return nothing
end   # function 
#----------------------------------------------------------------------
function write_params(fid::IOStream,params::Vector{Float64},pard::Vector{String},parc::Vector{String})

      hdr=" PARAMETERS FOLLOW";
      space5=blanks(5);             # arbitrary
      nparams=length(params);
      @printf(fid,"%s%5i%s\n",space5,nparams,hdr);

      space2=blanks(2);
      space4=blanks(4);
      for i=1:nparams
        pn = @sprintf("P%3.3i",i);
        l1=length(pard[i]);
        spacen=blanks(16-l1);
        description = @sprintf("%s%s%s%s%s%s",space2,pn,space2,pard[i],spacen,parc[i]);

        par=params[i];
        @printf(fid,"%s%14.6e%s\n",space2,par,description);
      end  

      return nothing
end   # function
#---------------------------------------------------------------------- 

function write_ps_hdr(fid::IOStream,nps::Int64)
#     This needs to be modified in a proper implementation

      hdr="Lines of passive scalar data follows: CONDUCT; RHOCP";
      space2=blanks(2);
      space5=blanks(5);
      @printf(fid,"%s%2i%s%s\n",space5,nps,space2,hdr);

      return nothing
end 
#----------------------------------------------------------------------
function write_logicalswitches_hdr(fid::IOStream,nlogic::Int64)

      # header
      space2=blanks(2);
      space5=blanks(5);
      hdr="LOGICAL SWITCHES FOLLOW";
      @printf(fid,"%s%3i%s%s\n",space5,nlogic,space2,hdr);

      return nothing
end   # function
#---------------------------------------------------------------------- 
function write_prenek_hdr(fid::IOStream)

      space5 = blanks(5);
      space4 = blanks(4);
      hdr    = "XFAC,YFAC,XZERO,YZERO";

      xfac  = 1.0
      yfac  = 1.0
      xzero = 0.0
      yzero = 0.0
      @printf(fid,"%10f%s%10f%s%10f%s%10f%s %s\n",xfac,space4,yfac,space4,xzero,space4,yzero,space5,hdr);
      
      return nothing
end
#---------------------------------------------------------------------- 
function write_mesh_hdr(fid,ndim,nelg,nelv,ifre2)

      if ndim==3
        hdr=" **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8"
        @printf(fid,"%s\n",hdr);
      else
        hdr=" **MESH DATA** 2 lines are X,Y. Columns corners 1-4";
        @printf(fid,"%s\n",hdr);
      end  

      space11=blanks(11);
      hdr="NEL,NDIM,NELV";
      if ifre2
        nelg=-nelg;
      end  
      @printf(fid,"%12i%3i%12i%s%s\n",nelg,ndim,nelv,space11,hdr); 

end
#----------------------------------------------------------------------
function write_mesh_coords(fid,ndim,XC,YC,ZC)

      @assert size(XC,2) == size(YC,2) "Unequal no of elements in XC and YC"

      if ndim == 3
        @assert size(XC,2) == size(ZC,2) "Unequal no of elements in XC and ZC"
      end  

      # @printf "Writing mesh"

      nelg=size(XC,2);

      # Mesh data   
      zlev = 1;      # Ignoring "zlevels"
      a = " ";
      igroup = 0;
      for e in 1:nelg 
        write_element_hdr(fid,e,zlev,a,igroup);
        if (ndim == 2)
          write_element_xyz(fid,XC[:,e],YC[:,e])
        else
          write_element_xyz(fid,XC[:,e],YC[:,e],ZC[:,e])
        end
      end

      return nothing
end   
#---------------------------------------------------------------------- 
function write_element_hdr(fid,e,zlev,a,igroup)

      space5 = blanks(5);
      space4 = blanks(4);
         
      @printf(fid,"%sELEMENT%12i [%5i%1s]%sGROUP%5i\n",space5,e,zlev,a,space4,igroup);

      return nothing
end   # function 
#----------------------------------------------------------------------  
function write_element_xyz(fid::IOStream,x::Vector{Float64},y::Vector{Float64})

      f = Ref(Printf.Format("%14.6E"));
      ind = 1:4
      Printf.format.(fid,f,x[ind])
      @printf(fid,"\n")
      Printf.format.(fid,f,y[ind])
      @printf(fid,"\n")
      
      return nothing
end
#----------------------------------------------------------------------
function write_element_xyz(fid::IOStream,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})

      f = Ref(Printf.Format("%14.6E"));
      ind = 1:4
      Printf.format.(fid,f,x[ind])
      @printf(fid,"\n")
      Printf.format.(fid,f,y[ind])
      @printf(fid,"\n")
      Printf.format.(fid,f,z[ind])
      @printf(fid,"\n")

      ind = 5:8
      Printf.format.(fid,f,x[ind])
      @printf(fid,"\n")
      Printf.format.(fid,f,y[ind])
      @printf(fid,"\n")
      Printf.format.(fid,f,z[ind])
      @printf(fid,"\n")
      
      return nothing
end
#----------------------------------------------------------------------
function write_element_xyz(fid::IOStream,ndim::Int64,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})

      f = Ref(Printf.Format("%14.6E"));
      if ndim==2
        ind = 1:4
        Printf.format.(fid,f,x[ind])
        @printf(fid,"\n")
        Printf.format.(fid,f,y[ind])
        @printf(fid,"\n")
      else
        ind = 1:4
        Printf.format.(fid,f,x[ind])
        @printf(fid,"\n")
        Printf.format.(fid,f,y[ind])
        @printf(fid,"\n")
        Printf.format.(fid,f,z[ind])
        @printf(fid,"\n")

        ind = 5:8
        Printf.format.(fid,f,x[ind])
        @printf(fid,"\n")
        Printf.format.(fid,f,y[ind])
        @printf(fid,"\n")
        Printf.format.(fid,f,z[ind])
        @printf(fid,"\n")
      end  
      
      return nothing
end
#----------------------------------------------------------------------
function write_curve_hdr(fid,ncurve)

      hdr = blanks(1)*"***** CURVED SIDE DATA *****"
      @printf(fid,"%s\n",hdr);
      hdr = "Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE";
      space2=blanks(2);
      @printf(fid,"%12i%s%s\n",ncurve,space2,hdr);

      return nothing
end 
#----------------------------------------------------------------------
function write_curvedata(fid::IOStream,nelg::Int64,ncurve::Int64,curveieg::Vector{Int64},curveedge::Vector{Int64},curveparams::Array{T},curvetype::Vector{String}) where {T<:AbstractFloat}

      datafmt = "%14.6E"^5
      space1  = blanks(1)
      if nelg<1000
        str = "%3i%3i"*datafmt*space1*"%1s\n"
        fmt = Printf.Format(str)
      elseif nelg<10^6
        str = "%2i%6i"*datafmt*space1*"%1s\n"
        fmt = Printf.Format(str)
      else
        str = "%2i%12i"*datafmt*space1*"%1s\n"
        fmt = Printf.Format(str)
      end

      for e in 1:ncurve
        ieg  = curveieg[e]
        edge = curveedge[e]
        cp1  = curveparams[1,e]
        cp2  = curveparams[2,e]
        cp3  = curveparams[3,e]
        cp4  = curveparams[4,e]
        cp5  = curveparams[5,e]
        ct   = curvetype[e]

        Printf.format(fid,fmt,edge,ieg,cp1,cp2,cp3,cp4,cp5,ct)
      end        

#%     60 format(i3,i3,5g14.6,1x,a1)
#%     61 format(i2,i6,5g14.6,1x,a1)
#%     62 format(i2,i12,5g14.6,1x,a1)
      return nothing

end   
#---------------------------------------------------------------------- 

function write_fluidBC_hdr(fid)

      hdr = blanks(2)*"***** BOUNDARY CONDITIONS ****";
      @printf(fid,"%s\n",hdr);
      hdr = blanks(2)*"***** FLUID   BOUNDARY CONDITIONS *****";
      @printf(fid,"%s\n",hdr);

      return nothing
end
#---------------------------------------------------------------------- 
function write_fluidBC(fid,ndim,nelg,BC0,Par0)

#     Either there is a bug for large element numbers or 
#     there is an inherent ordering of elements for BCs.

      # datafmt  = repmat('%14.6f',1,5);
      # f1       = Ref(Printf.Format("%14.6f"));
      # datafmt2 = repmat('%18.11f',1,5);
      # f2       = Ref(Printf.Format("%18.11f"));

      bc   = "%3s"           # Boundary Condition
      p1   = "%14.6f"^5      # Parameters
      p2   = "%18.11f"^5     # Parameters
      bl   = blanks(1)       # Blank
      br   = "\n"            # Break
      fmt1 = false

      if nelg<1000      
        #fmt = [space1,'%3s%3i%3i',datafmt,'\n'];
        ef  = "%3i%3i"         # El, Face
        str = bl*bc*ef*p1*br
        fmt = Printf.Format(str);
        fmt1 = true;
      elseif nelg<10^5
        #fmt = [space1,'%3s%5i%1i',datafmt,'\n'];
        ef  = "%5i%1i"         # El, Face
        str = bl*bc*ef*p1*br
        fmt = Printf.Format(str);

        fmt1 = true;

      elseif nelg<10^6
        #fmt = [space1,'%3s%6i',datafmt,'\n'];
        ef  = "%6i"         # El
        str = bl*bc*ef*p1*br
        fmt = Printf.Format(str);

        fmt1 = false
      else
        #fmt = [space1,'%3s%12i',datafmt2,'\n'];
        ef  = "%12i"         # El
        str = bl*bc*ef*p2*br
        fmt = Printf.Format(str);
       
        fmt1 = false
      end

      nfaces=2*ndim;
      for e in 1:nelg
        for j in 1:nfaces
          BC=BC0[j,e];
          c2=Par0[1,j,e] # mesh.cbc(j,e).connectsto;
          of=Par0[2,j,e] # mesh.cbc(j,e).onface;
          p3=Par0[3,j,e] # Par 3
          p4=Par0[4,j,e] # Par 4
          p5=Par0[5,j,e] # Par 5

          if (fmt1)
            Printf.format(fid,fmt,BC,e,j,c2,of,p3,p4,p5);
          else
            Printf.format(fid,fmt,BC,e,c2,of,p3,p4,p5);
          end
        end
      end  

#     20 FORMAT(1x,A3,2I3,5G14.6)
#     21 FORMAT(1x,A3,i5,i1,5G14.6)
#     22 FORMAT(1x,A3,i6,5G14.6)
#     23 FORMAT(1x,A3,i12,5G18.11)

      return nothing
end
#----------------------------------------------------------------------
function write_zerothermalBC_hdr(fid::IOStream)

      hdr = blanks(1)*"***** NO THERMAL BOUNDARY CONDITIONS *****"
      @printf(fid,"%s\n",hdr);

      return nothing
end 
#----------------------------------------------------------------------
function write_restart(fid::IOStream,nrst::Int64,rstfiles::Vector{String},rstoptions::Vector{String})

      space5 = blanks(5)
      space1 = blanks(1)
      hdr="PRESOLVE/RESTART OPTIONS"

      @printf(fid,"%s%3i%s%s\n",space5,nrst,space1,hdr);

      for i in 1:nrst
        @printf(fid,"%s%s%s\n",rstfiles[i],space1,rstoptions[i]);
      end

      return nothing
end
#---------------------------------------------------------------------- 
function write_ic(fid::IOStream,nic::Int64,icfiles::Vector{String})

      space5 = blanks(5);
      space1 = blanks(1);
      hdr    = "INITIAL CONDITIONS"

      @printf(fid,"%s%3i%s%s\n",space5,nic,space1,hdr);

      for i in 1:nic
        @printf(fid,"%s\n",icfiles[i]);
      end  

      return nothing
end 
#---------------------------------------------------------------------- 
function write_driveforce(fid::IOStream,ndf::Int64,driveforce::Vector{String})

      hdr     = blanks(1)*"***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q"
      @printf(fid,"%s\n",hdr);

      space5 = blanks(5);
      space1 = blanks(1);
      hdr    = "Lines of Drive force data follow";

      @printf(fid,"%s%3i%s%s\n",space5,ndf,space1,hdr);

      for i in 1:ndf
        @printf(fid,"%s\n",driveforce[i]);
      end  
      
      return nothing
end
#---------------------------------------------------------------------- 
function write_varprop(fid::IOStream,nvp::Int64,npackets::Int64,datapackets::Vector{String})

      hdr= blanks(1)*"***** VARIABLE PROPERTY DATA ***** Overrides Parameter data "
      @printf(fid,"%s\n",hdr)

      space5 = blanks(5)
      space1 = blanks(1)

      hdr = "Lines follow"
      @printf(fid,"%s%3i%s%s\n",space5,nvp,space1,hdr)

      hdr = "PACKETS OF DATA FOLLOW"
      @printf(fid,"%s%3i%s%s\n",space5,npackets,space1,hdr)
     
      for i in 1:npackets
        @printf(fid,"%s\n",datapacket[i]);
      end  

      return nothing
end
#---------------------------------------------------------------------- 
function write_hist(fid::IOStream,nhist::Int64,history::Vector{String})

      hdr=blanks(1)*"***** HISTORY AND INTEGRAL DATA *****"
      @printf(fid,"%s\n",hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      hdr = "POINTS. Hcode, I,J,H,IEL";
      @printf(fid,"%s%3i%s%s\n",space5,nhist,space1,hdr);
     
      for i=1:nhist
        @printf(fid,"%s\n",history[i]);
      end  

      return nothing
end
#---------------------------------------------------------------------- 
function write_outputspec(fid::IOStream,niospec::Int64,flds::Vector{String},Spec::Vector{Integer})

      hdr = blanks(1)*"***** OUTPUT FIELD SPECIFICATION *****"
      @printf(fid,"%s\n",hdr)

      space5 = blanks(5)
      space1 = blanks(1)

      hdr = "SPECIFICATIONS FOLLOW"
      @printf(fid,"%s%3i%s%s\n",space5,niospec,space1,hdr)

      space2 = blanks(2)
      for i in 1:niospec
        @printf(fid,"%s%i%s%s\n",space2,Spec[i],space5,flds[i]);
      end  

      return nothing
end
#---------------------------------------------------------------------- 
function write_objectspec(fid::IOStream,nobj::Int64,Objects::Vector{String},ObjSpec::Vector{Int64})


      hdr = blanks(1)*"***** OBJECT SPECIFICATION *****"
      @printf(fid,"%s\n",hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      space2 = blanks(2);
      for i in 1:nobj
        @printf(fid,"%s%3i%s%s\n",space2,ObjSpec[i],space5,Objects[i]);
      end  

      return nothing
end
##---------------------------------------------------------------------- 
function write_objectspec(fid::IOStream,Obj::NekObjects)


      hdr = blanks(1)*"***** OBJECT SPECIFICATION *****"
      @printf(fid,"%s\n",hdr);

      space5 = blanks(5);
      space1 = blanks(1);

      space2 = blanks(2);
      @printf(fid,"%s%3i%s%s\n",space2,Obj.surface,space5,"Surface Objects"); 
      @printf(fid,"%s%3i%s%s\n",space2,Obj.volume,space5, "Volume  Objects");
      @printf(fid,"%s%3i%s%s\n",space2,Obj.edge,space5,   "Edge    Objects");
      @printf(fid,"%s%3i%s%s\n",space2,Obj.point,space5,  "Point   Objects");

      return nothing
end
#---------------------------------------------------------------------- 

#function write_logicalswitches(fid,rea)
#
#      nlogic=rea.Nlogical;
#
#      # header
#      space2=blanks(2);
#      space5=blanks(5);
#      hdr="LOGICAL SWITCHES FOLLOW";
#      @printf(fid,"%s%3i%s%s\n",space5,nlogic,space2,hdr);
#
#      # Switches      
#      for i=1:nlogic
#        flag=rea.logical[i,2];
#        val =rea.logical[i,1];            
#        l1=length(val);
#        if l1==1
#          @printf(fid,"%s%s%s%s\n",space2,val,space5,flag);
#        else
#          val2=[];
#          for j=1:l1
#            val2 = val2*val[j]*" ";
#          end
#          val3 = val2[1:end-1]
#          @printf(fid,"%s%s%s%s\n",space2,val3,space2,flag);
#        end  
#      end  
#
#      return nothing
#end   # function
#---------------------------------------------------------------------- 

function blanks(n::Int)

  b = ""
  for i in 1:n
    b = b*" "
  end    

  return b
end
#---------------------------------------------------------------------- 



