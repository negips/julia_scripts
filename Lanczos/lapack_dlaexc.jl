function dlaexc!(wantq::Int, j1::Int, n1::Int, n2::Int, T::AbstractMatrix{Float64}, Q::AbstractMatrix{Float64})

  LinearAlgebra.chkstride1(T, Q)
  n = LinearAlgebra.checksquare(T)
  ldt = max(1, stride(T, 2))
  ldq = max(1, stride(Q, 2))
  work = Vector{Float64}(undef, n)
  info = Ref{Int}()
  #info = Vector{Int}(undef,1)    
  ccall((:dlaexc_, localsrc), Cvoid,
        ( Ref{Int},Ref{Int},
          Ptr{Float64},Ref{Int},Ptr{Float64},Ref{Int},
          Ref{Int},Ref{Int},Ref{Int},
          Ptr{Float64},Ptr{Int}
         ),
         wantq,n,
         T,ldt,Q,ldq,
         j1,n1,n2,
         work,info
       )

  #@show info
  #LinearAlgebra.LAPACK.chklapackerror(info[])
  #T, Q
  #return info
end





