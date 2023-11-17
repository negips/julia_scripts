#!/usr/bin/julia
#println("Testing Gather Scatter library wrapper")

using MPI

nid0  = 0

if (MPI.Initialized() == false)
  MPI.Init()
end  

comm = MPI.Comm
world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)

world2 = MPI.Comm_dup(world)



#println("We are here in $rank")



