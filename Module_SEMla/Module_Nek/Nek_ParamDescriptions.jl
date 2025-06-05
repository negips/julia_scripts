function param_descriptions()

  p = Vector{String}(undef,200)
  c = Vector{String}(undef,200)

  # This has been relunctantly compiled
  p[1]   = "DENSITY"
  c[1]   = "(For constant property case)"

  p[2]   = "VISCOSITY"
  c[2]   = "(==Re if <0. Otherwise ==1/Re)"

  p[3]   = "BETAG"
  c[3]   = "(Not in use)"

  p[4]   = "GTHETA"
  c[4]   = "(Not in use)"

  p[5]   = "PGRADX"
  c[5]   = "(Not in use)"
  
  p[6]   = "-"
  c[6]   = "(Not in use)"

  p[7]   = "RHOCP"
  c[7]   = "(Not in use)"

  p[8]   = "CONDUCT"
  c[8]   = "(Conductivity for constant properties. ==Pe if <0)"

  p[9]   = "-"
  c[9]   = "(Not in use. CPFLD(2,3)=PARAM(9))"

  p[10]  = "FINTIME"
  c[10]  = "(Simulation end time if >0)"

  p[11]  = "NSTEP"
  c[11]  = "(Number of time steps)"

  p[12]  = "DT"
  c[12]  = "(if <0 dt=|PARAM(12)|. Else dt upper bound)"

  p[13]  = "IOCOMM"
  c[13]  = "(Frequency of communication histories)"

  p[14]  = "IOTIME"
  c[14]  = "(if >0 time interval to dump fld file)"

  p[15]  = "IOSTEP"
  c[15]  = "(No. of time steps between fld dumps)"

  p[16]  = "PSSOLVER"
  c[16]  = "(Heat/PS solver. 1: Helmholtz 2:CVODE 3:CVODE with user jacobian)"

  p[17]  = "AXIS"
  c[17]  = "(Not in use)"

  p[18]  = "GRID"
  c[18]  = "(Not in use)"

  p[19]  = "INTYPE"
  c[19]  = "(Not in use)"

  p[20]  = "NORDER"
  c[20]  = "(Not in use)"

  p[21]  = "DIVERGENCE"
  c[21]  = "(Tolerance for pressure solver)"

  p[22]  = "HELMHOLTZ"
  c[22]  = "(Tolerance for velocity solver. Relative if <0)"

  p[23]  = "NPSCAL"
  c[23]  = "(No. of passive scalars)"

  p[24]  = "TOLREL"
  c[24]  = "(Relative tol for passive scalars)"

  p[25]  = "TOLABS"
  c[25]  = "(Absolute tol for passive scalars)"

  p[26]  = "COURANT"
  c[26]  = "(Maximum CFL. Or no of RK4 substeps for OIFS)"

  p[27]  = "TORDER"
  c[27]  = "(NBDINP. Temporal integration order)"

  p[28]  = "NABMSH"
  c[28]  = "(Temporal order for mesh velocity)"

  p[29]  = "MHD_VISCOS"
  c[29]  = "(==Magnetic Re if <0. Otherwise ==Magnetic viscosity if >0)"

  p[30]  = "USERVP"
  c[30]  = "(0: Constant properties 1,2: user-defined properties)"

  p[31]  = "NPERT"
  c[31]  = "(No. of perturbation modes. if >0 then time advance base flow.)"

  p[32]  = "NBCRE2"
  c[32]  = "(if >0: No. of BCs in re2 file. 0: All)"

  for i in 33:35
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[36]  = "XMAGNET"
  c[36]  = "(Not in use)"

  p[37]  = "NGRIDS"
  c[37]  = "(Not in use)"

  p[38]  = "NORDER2"
  c[38]  = "(Not in use)"

  p[39]  = "NORDER3"
  c[39]  = "(Not in use)"

  p[40]  = "-"
  c[40]  = "(Not in use)"

  p[41]  = "SEMG"
  c[41]  = "(1: multiplicative SEMG. (Not in use))"

  p[42]  = "PRES. SOLVER"
  c[42]  = "(0:GMRES 1:PCG)"

  p[43]  = "MULTILEVEL"
  c[43]  = "(0:Additive multilevel 1:Original two level)"

  p[44]  = "PRECONDITIONER"
  c[44]  = "(0:E-based additive Schwarz for PnPn-2 1:A-based)"

  p[45]  = "FACTOR"
  c[45]  = "(Free-surface stability control. Default 1.)"

  p[46]  = "IFSETICS"
  c[46]  = "(if >0 no call to SETICS)"

  p[47]  = "VNU"
  c[47]  = "(Poisson ratio for mesh elasticity solve)"

  p[48]  = "-"
  c[48]  = "(Not in use)"

  p[49]  = "TLFAC"
  c[49]  = "(Mixing length factor (Not in use))"

  p[50]  = "-"
  c[50]  = "(Not in use)"

  p[51]  = "-"
  c[51]  = "(Not in use)"

  p[52]  = "HISTEP"
  c[52]  = "(if >0 history points dump interval)"

  p[53]  = "-"
  c[53]  = "(Not in use)"

  p[54]  = "MFLOW"
  c[54]  = "(Direction of mass flow rate. 1:x, 2:y, 3:z)"

  p[55]  = "FLOWRATE"
  c[55]  = "(Volume flow rate for periodic case)"

  for i in 56:58
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[59]  = "EL. DEFORMED"
  c[59]  = "(if deformed elements, full jacobian evaluation)"

  p[60]  = "INITDISTB"
  c[60]  = "(Initialize velocity to 1e-10)"

  p[61]  = "-"
  c[61]  = "(Not in use)"
  
  p[62]  = "IFSWAP"
  c[62]  = "(Swap bytes for output)"

  p[63]  = "WDSIZO"
  c[63]  = "(output byte size. if >0 wdsizo=8)"

  p[64]  = "PERTRESTART"
  c[64]  = "(if =1 restart perturbation solution)"

  p[65]  = "IONODES"
  c[65]  = "(No. of IO nodes)"

  p[66]  = "OUT FORMAT"
  c[66]  = "(Output format. if >=0, binary. Else ASCII. Default: 6)"

  p[67]  = "RD FORMAT"
  c[67]  = "(Read format. if >=0, binary. Else ASCII)"

  p[68]  = "IASTEP"
  c[68]  = "(Averaging frequency in avg_all)"

  for i in 69:73
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[74]  = "HLMPRINT"
  c[74]  = "(If >0 print Helmholtz iterations)"

  for i in 75:81
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[82]  = "-"
  c[82]  = "(Coarse grid dimension (Not in use))"

  p[83]  = "-"
  c[83]  = "(Not in use)"

  p[84]  = "DT INIT"
  c[84]  = "(Force initial time step to this value)"

  p[85]  = "DT RATIO"
  c[85]  = "(dt=dtopf*PARAM(85))"

  p[86]  = "RESERVED"
  c[86]  = "(if >0 use skew-symmetric form for convection)"

  p[87]  = "-"
  c[87]  = "(Not in use)"

  p[88]  = "-"
  c[88]  = "(Not in use)"

  p[89]  = "RESERVED"
  c[89]  = "(Not in use)"

  p[90]  = "-"
  c[90]  = "(Not in use)"

  p[91]  = "-"
  c[91]  = "(Not in use)"

  p[92]  = "-"
  c[92]  = "(Not in use)"

  p[93]  = "MXPREV"
  c[93]  = "(No. of solutions for projection)"

  p[94]  = "MSTEP"
  c[94]  = "(Start of velocity projection)"

  p[95]  = "ISTART"
  c[95]  = "(Start of pressure projection)"

  p[96]  = "-"
  c[96]  = "(Not in use)"

  p[97]  = "-"
  c[97]  = "(Not in use)"

  p[98]  = "-"
  c[98]  = "(Not in use)"

  p[99]  = "DEALIASING"
  c[99]  = "(if <0: Disabled 3: Old dealiasing 4: New Dealiasing)"

  p[100] = "RESERVED"
  c[100] = "(Pres. preconditioner when using CG)"

  p[101] = "NMODES"
  c[101] = "(No. of additional modes for filtering)"

  p[102] = "-"
  c[102] = "(Not in use)"

  p[103] = "FIL. WGHT"
  c[103] = "(Filter weight)"

  for i in 104:106
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[107] = "ADDH2"
  c[107] = "(Add it to h2 in sethlm)"

  for i in 108:115
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  p[116] = "NELX"
  c[116] = "(No. of elements in x for Fast Tensor Product (FTP))"

  p[117] = "NELY"
  c[117] = "(No. of elements in y for Fast Tensor Product (FTP))"

  p[118] = "NELZ"
  c[118] = "(No. of elements in z for Fast Tensor Product (FTP))"

  for i=119:200
    p[i] = "-"
    c[i] = "(Not in use)"
  end  

  return p,c
end   # function
#---------------------------------------------------------------------- 

