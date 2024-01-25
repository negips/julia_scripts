# Lanczos Update

@doc """

        function BiOrthogonalize!(θl::AbstractVector{T},θr::AbstractVector{T},w::AbstractVector{T},v::AbstractVector{T},
                                  Ql::AbstractMatrix{T},Qr::AbstractMatrix{T},Rl::AbstractMatrix{T},Rr::AbstractMatrix{T},
                                  ngs::Int) where {T}

            Qr        - Orthonormal basis for Right Krylov Space
            Ql        - Orthonormal basis for Left Krylov Space
            Rr        - Triangular matrix for V = Qr*Rr
            Rl        - Triangular matrix for W = Ql*Rl
            v         - New right vector
            w         - New left vector
            θr        - Oblique Projections onto V along W
            θl        - Oblique Projections onto W along V
            ngs       - No of oblique Gram-Schmid procedures


"""
function BiOrthogonalize!(θl::AbstractVector{T},θr::AbstractVector{T},w::AbstractVector{T},v::AbstractVector{T},
                          Ql::AbstractMatrix{T},Qr::AbstractMatrix{T},Rl::AbstractMatrix{T},Rr::AbstractMatrix{T},
                          ngs::Int) where {T}


  @assert size(Ql,2) == size(Qr,2) "Unequal subspace sizes"

  zro  = T(0.0)
  m    = length(θl) 

  fill!(θl,zro)
  fill!(θr,zro)

  θv   = fill(zro,m)
  θw   = fill(zro,m)

  # New Lanczos Vectors
  for g = 1:ngs
    # Oblique Projections
    θw      = Rr'*Qr'*w
    θv      = Rl'*Ql'*v

    # Remove projections
    w       .= w  .- Ql*Rl*θw
    v       .= v  .- Qr*Rr*θv

    # Update total projections
    θl      .+= θw
    θr      .+= θv

  end     # g =1:ngs 

  return nothing
end

#----------------------------------------------------------------------
@doc """

      function BiorthoScale!(v::AbstractVector,w::AbstractVector)
      
      Scale w,v s.t <w,v> = 1.0.

"""
function BiOrthoScale!(w::AbstractVector{T},v::AbstractVector{T},x::AbstractVector{T}) where {T}

  w1      = copy(w)
  v1      = copy(v)

# Symmetric scaling
  θ       = w'*v                   # <w,v>
  δ       = sqrt(abs(θ))
  β       = (θ/δ)'

  v      .= v./δ
  w      .= w./β   

  αl      = (w'*w1)/(w'*w)
  αr      = (v'*v1)/(v'*v)

  return αl,αr
end
#---------------------------------------------------------------------- 

function GetQRviews(Q::AbstractMatrix,R::AbstractMatrix,k::Int)

  @assert k>=0
  
  Qview = view(Q,:,1:k)
  Rview = view(R,1:k,1:k)
  
  return Qview,Rview
end
#---------------------------------------------------------------------- 
@doc """

        function UpdateQR!(Q::AbstractMatrix{T},R::AbstractMatrix{T},x::AbstractVector{T},k::Int,ngs::Int) where {T}

            Q         - Orthonormal basis for the Krylov Space
            R         - Triangular matrix for V = Q*R
            x         - New right vector
            k         - Current size of Krylov space
            ngs       - No of oblique Gram-Schmid procedures

"""
function UpdateQR!(Q::AbstractMatrix{T},R::AbstractMatrix{T},x::AbstractVector{T},k::Int,ngs::Int) where {T}


  n         = length(x)
  θ         = zeros(T,k)
  QV        = view(Q,:,1:k)
  y         = view(Q,:,k+1)

  copyto!(y,1,x,1,n)

  for g = 1:ngs
    # Orthogonal Projections
    α        = QV'*y

    # Remove projections
    y       .= y  .- QV*α

    # Update total projections
    θ       .+= α
  end     # g =1:ngs

  β  = norm(y)
  y .= y/β
  r         = view(R,:,k+1)
  copyto!(r,1,θ,1,k)
  r[k+1] = β

  return nothing
end

#----------------------------------------------------------------------
@doc """

      SelectEigenvalues(evals::AbstractVector{T},Nev::Int) where {T}

      Output Int vector for selected eigenvalues.

"""
function SelectEigenvalues(evals::AbstractVector{T},Nev::Int) where {T}

      n           = length(evals)
      select      = fill(0,n)
      er          = real.(evals)
      ei          = imag.(evals)
      sind        = sortperm(ei,rev=true)   # Decreasing order of imaginary (eigenvalues).
      for i in 1:Nev
        j   = sind[i]
        select[j] = 1
      end  

      return select
end

#----------------------------------------------------------------------







