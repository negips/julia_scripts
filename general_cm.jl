#!/bin/julia

println("Testing Asymptotic Center-Manifold Matrix")

using LinearAlgebra
using SymPy

include("$JULIAHOME/Module_CenterManifold/CenterManifold.jl")
include("$JULIACOMMON/KroneckerDelta.jl")

l = 1
n = 1
m = 1
N = n+m
M = l+n+m

Ord = 2

# Symbols for all modes
Γstr        = "Γ"
λstr        = "λ"
Kstr        = "K"

#Cstr        = "ϕ"       # Critical modes
#Sstr        = "v"       # Stable modes
#Pstr        = "μ"       # Parameters
#Xstr        = "x"       # All modes
#Ystr        = "y"       # All Transformed modes
#θstr        = "θ"       # All eigenvalues
#QXstr       = "Q"       # All non-linearities in X
#QYstr       = "P"       # All non-linearities in Y
#αstr        = "α"       # Stable mode Approx. Coeffs
#
#C           = CenterManifold.GenerateSymbols(Cstr,NC)     # Critical Modes
#S           = CenterManifold.GenerateSymbols(Sstr,NS)     # Stable Modes
#μ           = CenterManifold.GenerateSymbols(Pstr,NP)     # Parameters
#X           = CenterManifold.GenerateSymbols(Xstr,nmodes) # All modes: NC + NP + NS. I just label them as "x"
#Y           = CenterManifold.GenerateSymbols(Ystr,nmodes) # All modes in the transformed system.
#θ           = CenterManifold.GenerateSymbols(θstr,nmodes) # Eigenvalue symbols for all modes

λ           = CenterManifold.GenerateSymbols(λstr,m)  # Oscillatory Eigenvalue symbols for all modes
Γ           = CenterManifold.GenerateSymMatrix(Γstr,n,l;ifreal=false)
Khat        = CenterManifold.GenerateSymMatrix(Kstr,M,M;ifreal=false)

#Khat        = Matrix{Sym}(undef,M,M)
#fill!(Khat,0)
#
#for i in 1:n
#  for j in 1:l
#    Khat[l+i,j] = Γ[i,j]
#  end
#end
#
#for i in 1:m
#  j = l+n+i
#  Khat[j,j] = λ[i]
#end

# Upper triangular elements = 0
for i in 1:M
  for j in i+1:M
    Khat[i,j] = 0
  end
end  

# First block = 0
for i in 1:l
  for j in 1:i-1
    Khat[i,j] = 0
  end
end  

# Oscillatory couling blocks = 0
for i in N+1:M
  for j in 1:i-1
    Khat[i,j] = 0
  end
end  

nt  = CenterManifold.NInteractionTerms(Ord,M)

Mat = fill(Sym(0),nt,nt)
# Second-Order Matrix
for i in 1:nt
  # x*x(i)
  ind_xx = CenterManifold.GetInteractionTerms(i,Ord,M)
  # display(ind_xx)
  for j in 1:nt
    ind_yy = CenterManifold.GetInteractionTerms(j,Ord,M)
    a = ind_yy[1]
    b = ind_yy[2]
    for k in 0:M-1
      ind1 = sort([k,b])
      if ind1 == ind_xx
        Mat[i,j] = Mat[i,j] + Khat[a+1,k+1]
      end  

      ind2 = sort([a,k])
      if ind2 == ind_xx
        Mat[i,j] = Mat[i,j] + Khat[b+1,k+1]
      end  
    end     # k
  end       # j
end         # i


# Mat4 = fill(Sym(0),nt,nt)
# # Second-Order Matrix
# for α in 0:M-1
#   for β in α:M-1
#     ind_xx = [α; β]
#     t      = CenterManifold.GetInteractionIndex(ind_xx,M)
# 
#     for i in 0:β
#       ind_yy    = [i; β]
#       k         = CenterManifold.GetInteractionIndex(ind_yy,M)
#       Mat4[t,k] = Mat4[t,k] + Khat[i+1,α+1]
#     end
# 
#     if (α != β)
#       for i in 0:α
#         ind_yy    = [i; α]
#         k         = CenterManifold.GetInteractionIndex(ind_yy,M)
#         Mat4[t,k] = Mat4[t,k] + Khat[i+1,β+1]
#       end
#     end  
# 
#     ##
#     ##  
#     for j in β:M-1
#       ind_yy    = [β; j]
#       k         = CenterManifold.GetInteractionIndex(ind_yy,M)
#       Mat4[t,k] = Mat4[t,k] + Khat[j+1,α+1]
#     end
# 
#     if (α != β)
#       for j in α:M-1
#         ind_yy    = [α; j]
#         k         = CenterManifold.GetInteractionIndex(ind_yy,M)
#         Mat4[t,k] = Mat4[t,k] + Khat[j+1,β+1]
#       end
#     end  
# 
#   end       # β 
# end         # α



Mat5 = fill(Sym(0),nt,nt)
# Second-Order Matrix
for α in 0:M-1
  for β in α:M-1
    ind_xx = [α; β]
    t      = CenterManifold.GetInteractionIndex(ind_xx,M)

    for i in 0:β
      ind_yy    = [i; β]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] + Khat[i+1,α+1]
    end

    for i in 0:α
      ind_yy    = [i; α]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] + Khat[i+1,β+1]
    end

    # Subtract α == β terms
    for i in 0:α
      ind_yy    = [i; α]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] - Khat[i+1,β+1]*KroneckerDelta(α,β)
    end
   
    ##
    ##  
    for j in β:M-1
      ind_yy    = [β; j]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] + Khat[j+1,α+1]
    end

    for j in α:M-1
      ind_yy    = [α; j]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] + Khat[j+1,β+1]
    end

    # Subtract α == β
    for j in α:M-1
      ind_yy    = [α; j]
      k         = CenterManifold.GetInteractionIndex(ind_yy,M)
      Mat5[t,k] = Mat5[t,k] - Khat[j+1,β+1]*KroneckerDelta(α,β)
    end

  end       # β 
end         # α

# Third Order
#---------------------------------------------------------------------- 

# Brute force
nt3    = CenterManifold.NInteractionTerms(3,M)
Mat_O3 = fill(Sym(0),nt3,nt3)
# Third Order Matrix
for j in 0:M-1
  for a in 0:M-1
    for b in a:M-1
      for c in b:M-1
        ind_yyy = [a;b;c]
        k       = CenterManifold.GetInteractionIndex(ind_yyy,M)

        ind_xxx = [j;b;c]
        t       = CenterManifold.GetInteractionIndex(ind_xxx,M)
        Mat_O3[t,k] = Mat_O3[t,k] + Khat[a+1,j+1]
       
        ind_xxx = [a;j;c]
        t       = CenterManifold.GetInteractionIndex(ind_xxx,M)
        Mat_O3[t,k] = Mat_O3[t,k] + Khat[b+1,j+1]

        ind_xxx = [a;b;j]
        t       = CenterManifold.GetInteractionIndex(ind_xxx,M)
        Mat_O3[t,k] = Mat_O3[t,k] + Khat[c+1,j+1]
      end   # c 
    end     # b
  end       # a
end         # j









