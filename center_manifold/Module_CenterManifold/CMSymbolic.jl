# I should turn this into a module when i'm done editing

#----------------------------------------------------------------------
""" 
      CollectTerms!(Collection::AbstractMatrix{T},G::AbstractVector{T},XNL::AbstractArray{T},d0::Dict) where {T}

      Collect terms of order N

"""
function CollectTerms!(Collection::AbstractArray{T},G::Union{AbstractVector{T},T},XNL::AbstractVector{T},d0::Dict) where {T}

  ns        = length(G)
  mt        = length(XNL)
  @assert ndims(XNL)==1 "XNL needs to be a Vector or a slice."

  s1        = Sym(1)

  for k in 1:ns
    Collection[:,k] = CollectTerms(G[k],XNL,d0)
  end

  return nothing
end
#---------------------------------------------------------------------- 
""" 
      CollectTerms(G::T,XNL::AbstractArray{T},d0::Dict) where {T}

      Collect terms of order N

"""
function CollectTerms(G::T,XNL::AbstractVector{T},d0::Dict) where {T}

  mt        = length(XNL)

  Collection = Vector{T}(undef,mt)

  for j in 1:mt
    a     = collect(G,XNL[j],evaluate=false)
    for (key,val) in a
      if key==XNL[j]
        Collection[j] = subs(val,d0)
        break
      else
        Collection[j] = 0
      end
    end
  end

  return Collection
end
#---------------------------------------------------------------------- 

""" 
      CollectOrderN!(Collection::AbstractMatrix{Sym},G::AbstractVector{Sym},XNL::Matrix{Sym},dict0::Vector{Dict},N::Int)

      Collect terms of order N

"""
function CollectOrderN!(Collection::AbstractArray{T},G::Union{AbstractVector{T},T},XNL::Matrix{T},dict0::Vector{Dict},N::Int) where {T}

  ns        = length(G)
  mt,mo     = size(XNL)
  dl        = length(dict0)
  @assert dl>=mo "Insufficient size of Dictionary dict0"
 
  for k in 1:ns
    fh        = expand(G[k])
    f1        = subs(fh,dict0[1])     # Get inhomogeneous part
    fh        = fh - f1               # fh is the homogeneous part
    for i in mo:-1:(N+1)
      f1      = subs(fh,dict0[i])     # Set ith order terms to zero
      fh      = copy(f1)
    end
    f1        = subs(fh,dict0[N])     # Set Nth order terms to zero
    fh        = fh - f1               # This gives us the Nth order terms

    nn        = length(dict0[N])
    for i in 1:nn
      # Replace the terms we want with 1.       
      # So they turn into inhomogeneous terms. Expensive procesure. 
      f1              = subs(fh,XNL[i,N]=>1)

      # Zero all homogeneous terms.
      Collection[i,k]   = subs(f1,dict0[1])        
    end           # i
  end             # k

  return nothing
end
#---------------------------------------------------------------------- 
@doc """
      GenerateMR(Collection::AbstractArray{T},α::AbstractArray{T}) where {T}

      Generate Matrix and RHS after Collection of terms in the Center Manifold Approximation.
      By Creating a Matrix we implcitly assume that we have linear terms

"""
function GenerateMR(Collection::AbstractArray{T},α::AbstractArray{T}) where {T}

  s  = size(α)
  ns = s[1]
  s  = size(Collection)
  nc = s[1]

  local Mat = Matrix{T}(undef,nc,ns)
  local Rhs = Vector{T}(undef,nc)

  tupv = Vector{Tuple{T,Int}}(undef,ns)
  for i in 1:ns
    tupv[i] = tuple(α[i],0)
  end
  d0 = Dict(tupv)

  # Rhs
  for k in 1:nc
    fh = Collection[k]
    f1 = subs(fh,d0)
    Rhs[k] = -f1
  end

  v = Vector{Sym}(undef,ns)
  for k in 1:nc
    tr = Collection[k] + Rhs[k]
    CollectTerms!(v,tr,α,d0)
    for j in 1:ns
      Mat[k,j] = v[j]
    end
  end

  return Mat,Rhs
end
#---------------------------------------------------------------------- 
@doc """
      GetSolnExpressions(α,Eqs::AbstractVector{T},XCNL,NC::Int,NP::Int,NS::Int,ApxO::Int,dict1::Dict) where {T}

      Collect terms associated with each variable in XCNL from each of the equation expressions in Eqs.
      For each variable one obtains a matrix associated with the unknowns α and the right hand side.
      Due to the structure of Center-Manifold approximation that we have done, the matrix is always diagonal.
      Hence the α's can easily be computed (as expressions).

      The output is the Lhs (Diagonals), Rhs, and the solution α (as expressions).

"""
function GetSolnExpressions(α,Eqs::AbstractVector{T},XCNL::AbstractArray{T},NC::Int,NP::Int,NS::Int,ApxO::Int,dict1::Dict) where {T}

  CMD = NC+NP     # Center Manifold dimension

  zro = 0.0+im*0.0

  # α-solutions (as symbols)
  α_expr   = Array{T}(undef,size(α))
  fill!(α_expr,0)
  Rhs      = Matrix{Vector{T}}(undef,ApxO,NS)
  Lhs      = Matrix{Vector{T}}(undef,ApxO,NS)

  # Approximate "α" Order-by-Order
  for AOrd in 1:ApxO
  
    AOterms   = CenterManifold.NInteractionTerms(AOrd,CMD)
  
    if AOrd == 1
  
      # No LHS/RHS at first order
      for k in 1:NS
        Lhs[AOrd,k]    = [0]
        Rhs[AOrd,k]    = [0]
      end  
  
    else
      Xc        = Matrix{T}(undef,AOterms,NS)
      fill!(Xc,0)
      xnl       = view(XCNL,1:AOterms,AOrd)
      
      println("Getting terms of Order $(AOrd)")
      @views CenterManifold.CollectTerms!(Xc,Eqs,xnl,dict1)
        
      for k in 1:NS
        # Apparently the matrices are always diagonal.
        # This is the case only if the linear term is a diagonal matrix.
        println("Building Matrix and Rhs for Order $(AOrd) for k=$(k)")
        Mat,rhs = @views CenterManifold.GenerateMR(Xc[:,k],α[1:AOterms,AOrd,k])
        lhs     = copy(rhs)
        for i in 1:AOterms
          lhs[i] = Mat[i,i]
          # For now I'm not checking for zero lhs
          sol = rhs[i]/lhs[i]
          α_expr[i,AOrd,k] = sol
        end
        Lhs[AOrd,k] = lhs
        Rhs[AOrd,k] = rhs
      end     # k in 1:NS
    end       # AOrd == 1  
  end         # AOrd in 1:ApxO

  return Lhs,Rhs,α_expr
end  
#---------------------------------------------------------------------- 
@doc """
      GetCoeffs(Eqs::AbstractVector{T},X::AbstractVector{T}) where {T}

      Collect terms associated with each variable in XCNL from each of the equation expressions in Eqs.
      For each variable one obtains a matrix associated with the unknowns α and the right hand side.
      Due to the structure of Center-Manifold approximation that we have done, the matrix is always diagonal.
      Hence the α's can easily be computed (as expressions).

      The output is the Lhs (Diagonals), Rhs, and the solution α (as expressions).

"""
function GetCoeffs(Eq::T,X::AbstractVector{T}) where {T}

  nv = length(X)     # No of Variables

  zro = 0.0+im*0.0
  TS  = T("Temp")       # Temporary Symbol

  dT  = GenDict(X,TS)   # Map all variables to a temporary variable
  d0  = tuple(TS,0)

  EQ2 = copy(Eq)

  f   = subs(Eq,dT)     # Map everything to single variable
  f2  = subs(f,d0)      # Zero single variable. We are left with inhomogeneous parts
  EQ2 = EQ2 - f2        # Remove the inhomogeneous parts

  # Get highest order of terms
  if EQ2 != 0
    i    = 21     # assuming we don't go beyond 20
    EQ3  = subs(EQ2,dT)  # Map everything to single variable again
    dif  = 0*EQ3

    while (dif==0)
      i  = i-1
      TT = TS^i
      d  = tuple(TT,0)
      eq = subs(EQ3,d)
      dif = EQ3 - eq
    end
  else
    i = 0
  end

  if i==0
    return [],[],f2
  end

  XX,xx  = GetAllNonLinearTerms(X,i)
  β      = copy(XX)
  fill!(β,0)
  d0X = GenDict0(X)   # Map all X to 0

  for j in 1:i
    @views CollectTerms!(β[:,j],EQ2,XX[:,j],d0X)
  end  


  # Terms, Coefficients and inhomogeneous term
  return XX,β,f2
end  
#---------------------------------------------------------------------- 
@doc """
      GetPQMap(Eqs::AbstractVector{Sym},P::AbstractArray{Sym},X2Y::Dict{Sym,Sym},Y::AbstractVector{Sym},NLOrd::Int)

      Get Mapping of Non-Linear Coefficient Terms P->Q due to X->Y substitution

      Output: PQMap

"""
function GetPQMap(Eqs::AbstractVector{Sym},P::AbstractArray{Sym},X2Y::Dict{Sym,Sym},Y::AbstractVector{Sym},NLOrd::Int)

  nvars   = length(Y)
  neqs    = length(Eqs)

  PQMap   = Vector{Dict{Sym,Sym}}(undef,neqs)
  
  # Only the Non-Linear terms are Mapped.
  for k in 1:neqs
    mt         = CenterManifold.TotalInteractionTerms(NLOrd,nvars) - nvars
 
    # Substitute X -> Y
    tmpeq      = Eqs[k]
    tmpeq2     = expand(subs(tmpeq,X2Y))
 
    # Get Coefficients of all orders of Y
    pqv        = Vector{Tuple{Sym,Sym}}(undef,mt)
    y1,qq1,res = CenterManifold.GetCoeffs(tmpeq2,Y)
    r,c        = size(qq1)
    m          = 0
    for j in 2:c
      nt2 = CenterManifold.NInteractionTerms(j,nvars)
      for l in 1:nt2
        m      = m + 1
        pqv[m] = tuple(P[l,j,k],qq1[l,j])
      end
    end
    PQMap[k] = Dict(pqv)
  end

  return PQMap
end 
#---------------------------------------------------------------------- 
@doc """
      BuildGenericEquation(Qstr::String,NLO::Int,NEQ::Int,X::AbstractVector{T};nonlinearonly=false) where {T}

      Build a generic ̇∂X/∂t = Qi*Xi + Qij*Xi*Xj + Qijk*Xi*Xj*Xk ...

      with maximum order of terms = NLO.

      NEQ - No of Equations to construct

"""
function BuildGenericEquation(Qstr::String,NLO::Int,NEQ::Int,X::AbstractVector{T};nonlinearonly=false) where {T}

  # No of modes
  nm  = length(X)

  # Total interaction terms (of highest order)
  nt = NInteractionTerms(NLO,nm)

  # Create all Non-linear combination of terms
  XNL,xnl = GetAllNonLinearTerms(X,NLO)

  # Build Generic Coefficients
  Q  = Array{T}(undef,nt,NLO,NEQ)
  fill!(Q,0)

  for k in 1:NEQ
    if NEQ==1
      qq,q = GenerateAllNonLinearSymbols(Qstr,nm,NLO)
    else
      qq,q = GenerateAllNonLinearSymbols(Qstr,nm,NLO,super=k)
    end  
  
    Qvw = view(Q,:,:,k)
    copyto!(Qvw,1,qq,1,nt*NLO)
  end

  # Evolution equations 
  Xdot       = Vector{T}(undef,NEQ)
  fill!(Xdot,0)
  
  if nonlinearonly
    j1 = 2
  else
    j1 = 1
  end

  # Generic Non-Linear Evolution Equation
  for k in 1:NEQ
    for j in j1:NLO
      nt = NInteractionTerms(j,nm)
      for i in 1:nt
        Xdot[k] = Xdot[k] + Q[i,j,k]*XNL[i,j]
      end
    end  
  end


  return Q,XNL,Xdot
end  
#---------------------------------------------------------------------- 
@doc """
      AddLinearOP!(Xdot::AbstractVector{T},X::AbstractVector{T},LOP::AbstractMatrix{T}) where {T <: Number}

      Add a Linear Operator to the evolution ∂X/∂t = LOP*X + ∂X/∂t

"""
function AddLinearOP!(Xdot::AbstractVector{T},X::AbstractVector{T},LOP::AbstractMatrix{S}) where {T <: Number, S<:Number}

  # No of modes
  nm = length(X)
  ne = length(Xdot)

  r,c = size(LOP)

  @assert r == ne "length(Xdot) ≠ No. of rows(LOP)"
  @assert c == nm "length(X) ≠ No. of columns(LOP)"

  # Add a linear operator
  for k in 1:ne
    for j in 1:nm
      Xdot[k] = Xdot[k] + LOP[k,j]*X[j]
    end  
  end

  return nothing
end  
#---------------------------------------------------------------------- 

@doc """
      GetAllNonLinearTerms(X::AbstractVector{T},maxp::Int) where {T}

      Generate all non-linear combinations of terms in the vector X of 
      total order p = 1...maxp

      Terms output as a matrix.

"""
function GetAllNonLinearTerms(X::AbstractVector{T},maxp::Int) where {T}

  nm    = length(X)
  mmax  = NInteractionTerms(maxp,nm)
  mt    = TotalInteractionTerms(maxp,nm)

  XNL   = Matrix{T}(undef,mmax,maxp)
  fill!(XNL,0)
  xnl   = Vector{T}(undef,mt)
  fill!(xnl,0)

  i0 = 0
  for i in 1:maxp
    nt           = NInteractionTerms(i,nm)
    xx           = GetNonLinearTerms(X,i)
    XNL[1:nt,i]  = copy(xx)
    copyto!(xnl,i0+1,xx,1,nt)
    i0           = i0+nt
  end

  return XNL,xnl
end
#---------------------------------------------------------------------- 
@doc """
      GetNonLinearTerms(X::AbstractVector{T},maxp::Int) where {T}

      Generate all non-linear combinations of terms in the vector X for a total order
      of maxp.

"""
function GetNonLinearTerms(X::AbstractVector{T},maxp::Int) where {T}

  ns        = length(X)                    # No of symbols
  nt        = NInteractionTerms(maxp,ns)   # No of symbols
  NL        = Vector{T}(undef,nt) 

  p         = fill(0,1,ns)
  p[1,1]    = maxp
  i         = 1
  jlast     = 1

  # Put symbols in order
  nl        = InteractionNTuples(p,maxp,ns)
  NL[1]     = X[nl[1]]
  for i in 2:maxp
    NL[1] = NL[1]*X[nl[i]]
  end  
  
  while (i < nt) 
    i         = i + 1
    pw        = view(p,1,:)
    jlast = NextRow!(pw,maxp,jlast)

    # Put symbols in order
    nl        = InteractionNTuples(p,maxp,ns)
    NL[i]     = X[nl[1]]
    for j in 2:maxp
      NL[i] = NL[i]*X[nl[j]]
    end  

  end

  return NL
end
#----------------------------------------------------------------------
@doc """
      GenerateSymbols(sym::String,ns::Int;suffix="",prefix="",super::Int=-99,ifreal=false)

      Generate ns symbols with the string "sym". Numbering done in the subscript.
      Optionally adds a superscript string if provided.

"""
function GenerateSymbols(sym::String,ns::Int;suffix="",prefix="",super::Int=-99,ifreal=false)

  # ns       # No of symbols for "sym"
  SF        = SymbolFunction()

  S         = Vector{Sym}(undef,ns)
  if super!=-99
    str = sym*"^$(super)"
  else
    str = sym
  end

  if ns==1
    # No subscript if ns == 1
    str2    = prefix*str*suffix
    if (ifreal)
      S[1]  = SF(str2,real=true)
    else  
      S[1]  = SF(str2)
    end 
  else
    for i in 1:ns
      str2    = prefix*str*"_$(i)"*suffix
      if (ifreal)
        S[i]    = SF(str2,real=true)
      else  
        S[i]    = SF(str2)
      end  
    end
  end  

  return S
end
#---------------------------------------------------------------------- 
@doc """
      GenerateSymbolsSuper(sym::String,ns::Int;suffix="",prefix="",subs="",ifreal=false)

      Generate ns symbols with the string "sym". Numbering done in the superscript.
      Optionally adds a subscript string, suffix and prefix if provided.

"""
function GenerateSymbolsSuper(sym::String,ns::Int;suffix="",prefix="",subs="",ifreal=false)

  # ns       # No of symbols for "sym"
  SF        = SymbolFunction()
  S         = Vector{Sym}(undef,ns)
  if !isempty(subs)
    str = sym*"_"*subs*"^"
  else
    str = sym
  end  

  if ns==1
    # No superscript if ns == 1
    str2    = prefix*str*suffix
    if (ifreal)
      S[1]  = SF(str2,real=true)
    else  
      S[1]  = SF(str2)
    end 
  else
    for i in 1:ns
      str2    = prefix*str*"^$(i)"*suffix
      if (ifreal)
        S[i]    = SF(str2,real=true)
      else  
        S[i]    = SF(str2)
      end  
    end
  end  

  return S
end

#----------------------------------------------------------------------
@doc """
      GenerateSymMatrix(sym::String,ns::Int;super::Int=-99,ifreal=false)

      Generate Matrix of symbols with the string "sym". Numbering done in the subscript.
      Optionally adds a superscript string if provided.

"""
function GenerateSymMatrix(sym::String,nr::Int,nc::Int;super::Int=-99,ifreal=false)

  # ns       # No of symbols for "sym"
  SF        = SymbolFunction()

  S         = Matrix{Sym}(undef,nr,nc)
  if super!=-99
    str = sym*"^$(super)"
  else
    str = sym
  end

  for i in 1:nr
    for j in 1:nc
      str2     = str*"_$(i)_$(j)"
      if (ifreal)
        S[i,j] = SF(str2,real=true)
      else 
        S[i,j] = SF(str2)
      end  
    end
  end  

  return S
end
#---------------------------------------------------------------------- 

@doc """
      GenerateNonLinearSymbols(sym::String,ns::Int,maxp::Int;suffix="",prefix="",super::Int=-99,ifreal=false)

      Generate combinatorial symbols with the string "sym" for maximum power maxp.
      ns is the number of "linear" symbols. Numbering done in the subscript. 
      Optionally adds a superscript string if provided.

"""
function GenerateNonLinearSymbols(sym::String,ns::Int,maxp::Int;suffix="",prefix="",super::Int=-99,ifreal=false)

  #S   = GenerateSymbols(sym,ns)
  # NL  = GetNonLinearTerms(S,maxp)
  p   = fill(0,1,ns)
  nt  = NInteractionTerms(maxp,ns)

  SF  = SymbolFunction()
  NL  = Vector{Sym}(undef,nt)

  j   = 0
  if super!=-99
    str = sym*"^$(super)"
  else
    str = sym
  end  
 
  for i in 1:nt
    pw  = view(p,1,:)
    j   = NextRow!(pw,maxp,j)
    ind = InteractionNTuples(p,maxp,ns)
    strk = str
    for k in ind
      strk = strk*"_$k"
    end
    str2  = prefix*strk*suffix
    if (ifreal)
      NL[i]    = SF(str2,real=true)
    else  
      NL[i]    = SF(str2)
    end  
  end  

  return NL
end
#----------------------------------------------------------------------
@doc """
      GenerateAllNonLinearSymbols(sym::String,ns::Int,maxp::Int;super::Int=-99)

      Generate all combinatorial symbols with the string "sym" for total power p=1...maxp.
      ns is the number of "linear" symbols. Numbering done in the subscript. 
      Optionally adds a superscript string if provided.

"""
function GenerateAllNonLinearSymbols(sym::String,ns::Int,maxp::Int;super::Int=-99)

  nmax  = NInteractionTerms(maxp,ns)
  nt    = TotalInteractionTerms(maxp,ns)

  Q   = Array{Sym}(undef,nmax,maxp)
  fill!(Q,0)
  q   = Vector{Sym}(undef,nt)
  fill!(q,0)

  # Non-Linear Symbols
  i0 = 0
  for j in 1:maxp
    sup = super
    qnl = GenerateNonLinearSymbols(sym,ns,j;super=sup)
    l   = length(qnl)
    qvw = view(Q,:,j)
    copyto!(qvw,1,qnl,1,l)

    copyto!(q,i0+1,qnl,1,l)
    i0 = i0+l
  end  

  return Q,q
end


#----------------------------------------------------------------------
@doc """
      GenDict0(X::AbstractVector{T}) where {T}

      Generate Dictionary to map all symbols in X to zero

"""
function GenDict0(X::AbstractVector{T}) where {T}

   d0 = GenDict(X,0)

   return d0
end
#---------------------------------------------------------------------- 
@doc """
      GenDict(X::AbstractVector{T},V::S) where {T, S}

      Generate Dictionary to map all symbols in X to V

"""
function GenDict(X::AbstractVector{T},V::S) where {T, S}

   ns = length(X)
   tupv = Vector{Tuple{T,S}}(undef,ns)

   nonzero = 0
   for i in 1:ns
     s = X[i]
     if s != T(0)
       nonzero+=1
       tupv[nonzero] = tuple(s,V) 
     end  
   end

   return Dict(tupv[1:nonzero])
end
#---------------------------------------------------------------------- 
function SymbolFunction()

  # return Sym

  return SymPy.symbols

end
#---------------------------------------------------------------------- 



