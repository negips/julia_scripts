#!/usr/bin/julia
println("Testing Gather Scatter library wrapper")

using MPI

nid0  = 0

if (MPI.Initialized() == false)
  MPI.Init()
end  
  
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

N = 5;
gl_num = ones(Int64,N)

np = 1
gsh = 1

#ccall((:fgslib_gs_setup_, "./libgs.so"), Cvoid, (Int32, Ptr{Int64}, Int32, MPI.MPI_Comm, Int32), gsh, gl_num, N, comm, np)
comm = MPI.Comm
world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)

N = Int32(5);
glnum = zeros(Int64,N)

for i in 1:N-1
  glnum[i] = Int64(i)
end

glnum[N] = Int64(N-1)

gsh      = Int32(0)
np       = 1

#st = ccall((:fgslib_gs_setup_, "./libgs.so"), Cvoid, 
#           (Int32, Ptr{Int64}, comm, Int32), 
#           gsh, glnum, Ptr{world}, np)

function compute_dot(DX::Vector{Float64}, DY::Vector{Float64})
    @assert length(DX) == length(DY)
    n = length(DX)
    incx = incy = 1
    product = @ccall "liblapack".ddot_(
        n::Ref{Int32}, DX::Ptr{Float64}, incx::Ref{Int32}, DY::Ptr{Float64}, incy::Ref{Int32})::Float64
    return product
end

function cgs_setup(gsh::Int32,GLO_NUM::Vector{Int64},N::Int32,comm::MPI.Comm,NP::Int)
    @ccall "libgs".gslib_gs_setup(
       gsh::Ref{Int32}, GLO_NUM::Ptr{Int64}, N::Int32, comm::MPI.Comm, NP::Ref{Int64})::Nothing
end

# Fortran
function fgs_setup(gsh::Int32,GLO_NUM::Vector{Int64},N::Int32,comm,NP::Int)
    @ccall "./libgs.so".fgslib_gs_setup_(
      gsh::Ref{Int32},GLO_NUM::Ptr{Int64}, N::Ref{Int32}, comm::Ptr{MPI.MPI_Comm}, NP::Ref{Int64})::Nothing
end
#MPI.run_init_hooks()

#fgs_setup(gsh,glnum,N,world,np)







