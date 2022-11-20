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

ccall((:fgslib_gs_setup_, "./libgs.so"), Cvoid, (Int32, Ptr{Int64}, Int32, MPI.MPI_Comm, Int32), gsh, gl_num, N, comm, np)
comm = MPI.Comm
world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)

N = 5;
glnum = zeros(Int64,N)

for i in 1:N-1
  glnum[i] = Int64(i)
end

glnum[N] = Int64(N-1)

gsh      = 0
np       = 1

st = ccall((:fgslib_gs_setup_, "./libgs.so"), Cvoid, 
           (Int32, Ptr{Int64}, Int32, MPI.Comm, Int32), 
           gsh, glnum, N, Ptr{world}, np)



