#!/usr/bin/julia

println("Testing Julia interface for libgs")



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

