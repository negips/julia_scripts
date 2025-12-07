#
function Complex2String(val::ComplexF64)

  vr = real(val)
  vi = imag(val)

  str1      = ""
  space     = " "
  strr      = @sprintf("%.4f",abs(vr))
  stri      = @sprintf("%.4f",abs(vi))

  if (vr<0.0)
    if (vi<0.0)
      sni = "+"
    else
      sni = "-"
    end  
    str1 = str1*"- ("*strr*space*sni*space*stri*"im)"
  else
    if (vi<0.0)
      sni = "-"
    else
      sni = "+"
    end  
    str1 = str1*"+ ("*strr*space*sni*space*stri*"im)"
  end  

  return str1
end
#----------------------------------------------------------------------
function polyLatexString(var::String,pind::Vector{Int})

  strv = ""
  for k in 1:length(pind)
    strv = strv*var*"_{$(pind[k])}"
  end
  polystr = L"%$strv"

  return polystr
end  
#---------------------------------------------------------------------- 
function polystring(var::String,pind::Vector{Int})

  strv = ""
  for k in 1:length(pind)
    strv = strv*var*"$(pind[k])"
  end
  polystr = strv

  return polystr
end  
#---------------------------------------------------------------------- 

function CM_DisplayEquation(ind::Int,var::String,Khat::Matrix{T},G2::Matrix{T},G3::Matrix{T};tol=1.0e-10) where {T<:ComplexF64}

  nv,nv2    = size(Khat)

  str0      = "d"*var*"/dt"
  space     = " "

  ord = 1
  nt  = CenterManifold.NInteractionTerms(ord,nv)
  str1      = ""
  for i in 1:nt
    val     = Khat[ind,i]
    if (abs(val)>tol)
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str1    = str1*str*polystr*space
    end  
  end

  ord       = 2
  nt        = CenterManifold.NInteractionTerms(ord,nv)
  str2      = ""
  for i in 1:nt
    val     = G2[ind,i]
    if abs(val)>tol
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str2    = str2*str*polystr*space
    end  
  end

  ord       = 3
  nt        = CenterManifold.NInteractionTerms(ord,nv)
  str3      = ""
  for i in 1:nt
    val     = G3[ind,i]
    if abs(val)>tol
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str3    = str3*str*polystr*space
    end  
  end

  @printf("%s\n",str1)
  @printf("%s\n",str2)
  @printf("%s\n",str3)

  return nothing
end
#---------------------------------------------------------------------- 
function CM_DisplayTerms(ind::Int,var::String,Khat::Matrix{T},G2::Matrix{T},G3::Matrix{T};tol=1.0e-10) where {T<:ComplexF64}

  nv,nv2    = size(Khat)

  str0      = "d"*var*"/dt"
  space     = " "

  ord = 1
  nt  = CenterManifold.NInteractionTerms(ord,nv)
  for i in 1:nt
    val     = Khat[ind,i]
    if (abs(val)>tol)
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str1    = str*polystr
      @printf("%s\n",str1)
    end  
  end
  @printf("------------------------------\n")

  ord       = 2
  nt        = CenterManifold.NInteractionTerms(ord,nv)
  for i in 1:nt
    val     = G2[ind,i]
    if abs(val)>tol
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str2    = str*polystr
      @printf("%s\n",str2)
    end  
  end
  @printf("------------------------------\n")

  ord       = 3
  nt        = CenterManifold.NInteractionTerms(ord,nv)
  str3      = ""
  for i in 1:nt
    val     = G3[ind,i]
    if abs(val)>tol
      str     = Complex2String(val)
       
      pind    = CenterManifold.GetPolynomialIndices(i,ord,nv)
      pind    = pind .+ 1
      polystr = polystring(var,pind)
      str3    = str*polystr
      @printf("%s\n",str3)
    end  
  end
  @printf("------------------------------\n")

  return nothing
end
#---------------------------------------------------------------------- 






