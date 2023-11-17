#!/bin/julia
#println("Testing MPI in Julia")
using MPI
using ParallelUtilities

if !MPI.Initialized()
  MPI.Init()
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

# print("Hello world, I am rank $(rank) of $(size)\n")

test = 7 

# Scatter!
#--------------------------------------------------

if test == 1
  send_buf = nothing
  recv_buf = Vector{Float64}(undef,size)
  
  if rank == 0
    send_buf = rand(Float64,size,size)
    print("Scatter!: Original array on rank 0:\n $(send_buf)\n")
  end
  #rec_buf = MPI.Scatter(send_buf, Int, comm; root=0)
  MPI.Scatter!(send_buf, recv_buf, comm; root=0)
  
  print("I got this on Rank $(rank):\n $(recv_buf)\n")
end

# Scatterv!
#--------------------------------------------------  
if test == 2

  lengths  = rand(1:4,size)
  total    = sum(lengths)
  offsets  = cumsum(lengths) .- lengths
 
  MPI.Bcast!(lengths,comm; root=0)
  print("I got Lengths: $(lengths) on Rank $(rank):\n")

#  MPI.Bcast!(offsets,comm; root=0)
#  print("I got Offsets: $(offsets) on Rank $(rank):\n")

  send_buf = nothing
  recv_buf = Vector{Float64}(undef,lengths[rank+1])
 
  if rank == 0
    send_buf = rand(Float64,total)
    print("Scatterv!: Original array on rank 0:\n $(send_buf)\n")
  end
  #rec_buf = MPI.Scatter(send_buf, Int, comm; root=0)
  MPI.Scatterv!(VBuffer(send_buf, lengths, offsets, MPI.DOUBLE), recv_buf, comm; root=0)
  
  print("I got this on Rank $(rank):\n $(recv_buf)\n")
end



# Gather
#--------------------------------------------------

if test == 3
  send_buf = nothing
  v        = Vector{Float64}(undef,size)
 
  if rank == 0
    send_buf = rand(Float64,size,size)
    print("Scatter!: Original array on rank 0:\n $(send_buf)\n")
  end
  #rec_buf = MPI.Scatter(send_buf, Int, comm; root=0)
  MPI.Scatter!(send_buf, v, comm; root=0)

  print("I got this on Rank $(rank):\n $(v)\n")

  v = v.^2

#  recv_buf = MPI.Gather(v,comm; root=0)        # Only root gets it
  recv_buf = MPI.Allgather(v,comm)              # Everyone gets it

  print("I got this on Rank $(rank):\n $(recv_buf)\n")

end


# bcast
#--------------------------------------------------  
if test == 4

  send_buf = nothing

  if rank == 0
    send_buf = "Hello from Rank 0."
  end

  recv_buf = MPI.bcast(send_buf, comm; root=0)
  print("Received Message: $(recv_buf) on Rank $(rank)\n")
 
end


# Reduce
#--------------------------------------------------  
if test == 5

  length = 5
  v = rand(Float64,length)
  print("Original vector on Rank $(rank):\n $(v)\n")

  recv_buf = MPI.Reduce(v, +, comm; root=0)
  
  if rank == 0
    print("New Array on rank 0\n $(recv_buf)\n")
    total = sum(recv_buf)
    print("Total sum: $(total)\n")
  end  

end

# Reduce with overloading +
#--------------------------------------------------
if test == 6

  import Base.+

  struct Point{T}
    x::T
    y::T
  end  

  +(A::Point{T}, B::Point{T}) where T = Point{T}(A.x + B.x, A.y + B.y)

  p = Point(rand(),rand())
#  p = Point(rank,rank)

  print("Original point on Rank $(rank):\n $(p)\n")

  recv_buf = MPI.Reduce(p, +, comm; root=0)
  
  if rank == 0
    print("New point on rank 0\n $(recv_buf)\n")
  end  

end



# Point to Point Communication
#-------------------------------------------------- 

# send/recv
#--------------------------------------------------
if test == 7

  import Base.+

  struct Point{T}
    x::T
    y::T
  end  

#  +(A::Point{T}, B::Point{T}) where T = Point{T}(A.x + B.x, A.y + B.y)

  p = Point(rand(),rand())
#  p = Point(rank,rank)

#  if mod(rank,2) == 0
#    destn = mod(rank+1,size)
#    if destn != 0
#      MPI.send(p,comm,dest=destn)
#      print("Original point data sent from Rank $(rank) to Rank $(destn):\n $(p)\n")
#    end
#  else
#    src = rank-1
#    data = MPI.recv(comm,source=src)
#    print("point data received on Rank $(rank):\n $(data)\n")
#  end

  comm2 = MPI.Comm_dup(comm)

  destn = (rank+1)%size
  MPI.send(p,comm2,dest=destn)
  print("Original point data sent from Rank $(rank) to Rank $(destn):\n $(p)\n")
  src = rank-1
  if src<0
    src = size-1
  end
  data = MPI.recv(comm2,source=src)
  print("point data received on Rank $(rank) from Rank $(src):\n $(data)\n")

end



# Isend/Irecv!
#--------------------------------------------------
if test == 8


  length   = 2
  send_buf = rand(Float64,length) 
  recv_buf = Vector{Float64}(undef,length)

  destn = (rank+1)%size
  send_status = MPI.Isend(send_buf,comm;dest=destn)

  src = rank-1
  if src<0
    src = size-1
  end
  recv_status = MPI.Irecv!(recv_buf,comm; source=src)

  stat =  MPI.Wait!(recv_status)

  print("Original data sent from Rank $(rank) to Rank $(destn):\n $(send_buf)\n")
  print("Data received on Rank $(rank) from Rank $(src) (Status=$(stat.error)):\n $(recv_buf)\n")

  MPI.Barrier(comm)
end



# Reduce with (random) view
#--------------------------------------------------  
if test == 9

  VT     = Float64
  length = 5
  vlen   = 3
  vw     = rand(1:length,vlen)
  A      = rand(VT,length)
  B      = Vector{VT}(undef,vlen)

  Avw    = view(A,vw)
  copyto!(B,1,Avw,1,vlen)

  print("Original vector on Rank $(rank):\n $(A)\n")

  recv_buf = MPI.Allreduce(B, +, comm)
 
  copyto!(Avw,1,recv_buf,1,vlen)

  print("New Array on rank $(rank):\n $(A)\n")

  if rank == 0
    total = sum(recv_buf)
    print("Total sum: $(total)\n")
  end  

end


# Generate a random list of numbers
#--------------------------------------------------  
if test == 10
  n = 5
  lst = rand(1:n, n)
  
  # Perform parallel list ranking
  ranks = parallel_list_ranking(lst)
  
  # Print the resulting ranks
  println(ranks)
end  



#---------------------------------------------------------------------- 
MPI.Finalize()











