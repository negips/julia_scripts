println("Testing mxm calls from Julia")

using Random
using LinearAlgebra, StaticArrays
using BenchmarkTools, TimerOutputs


function stupidf_julia(a::Float64,b::Float64)

      c = ccall((:stupidf_, "./mxm"),
                Float64,
                (Ref{Float64}, Ref{Float64}),
                a,b)


      return c
end

function stupids_julia!(a::Float64,b::Float64,c::Float64)


      c1 = Ref(c)

      ccall((:stupids_wrapper, "./mxm"),
                Cvoid,
                (Ref{Float64}, Ref{Float64}, Ref{Float64}),
                a,b,c1)

      println(a)
      println(b)
      println(c)
      
end

function mxm_func_julia(a,n1::Int64,b,n2::Int64,n3::Int64)
#     mxm(a,n1,b,n2,c,n3)
#     Compute matrix-matrix product C = A*B
#     for contiguously packed matrices A,B, and C.

      a1 = Ref(a)
      b1 = Ref(b)

      println([n1 n2 n3])

      c = ccall((:fmxm_, "./mxm"),Float64,
                (Ptr{Float64}, Ref{Int64}, Ptr{Float64}, Ref{Int64}, Ref{Int64}),
                a,n1,b,n2,n3)

      println(c)
      
      return c
end

function mxm_sub_julia!(a,n1,b,n2,c,n3)
#     mxm(a,n1,b,n2,c,n3)
#     Compute matrix-matrix product C = A*B
#     for contiguously packed matrices A,B, and C.

ccall((:mxm_, "./mxm"),Matrix{Float64},(Ref{Float64},Ref{Int64},Ref{Float64},Ref{Int64},Ref{Float64},Ref{Int64}),a,n1,b,n2,c,n3)

      return c
end


const to = TimerOutput()
const al = TimerOutput()

nr = 5
nc = 5

nr1 = 5

# A  = rand(Float64, (nr, nc))
# B  = rand(Float64, (nr, nc))
# C  = rand(Float64, (nr, nc))
# Wk = zeros(Float64, (nr, nc))

# A  = @timeit al "rand" @SMatrix rand(nr, nc)
# B  = @timeit al "rand" @SMatrix rand(nr, nc)
# C  = @timeit al "rand" @SMatrix rand(nr, nc)
# Wk = @SMatrix zeros(nr, nc)

niter = 1


for i in 1:niter

   A  = @timeit al "rand" rand(Float64,(nr, nc))
   B  = @timeit al "rand" rand(Float64,(nr, nc))
   C  = @timeit al "rand" rand(Float64,(nr, nc))
   D  = @timeit al "rand" rand(Float64,(nr, nc))
   E  = @timeit al "rand" rand(Float64,(nr, nc))

   global A, B, C, D, E 

#   C = A*B
   @timeit to "mult" mul!(C,A,B)
   C = @timeit to "std" A*B 
#   D = @timeit to "mxm" mxm_sub_julia!(A,nr,B,nr,D,nc)
   E = @timeit to "mxm" mxm_func_julia(A,nr1,B,nr,nc)

end

a = 1.0
b = 2.0
e = 0.0

d = @timeit to "stu" stupidf_julia(a,b)
@timeit to "stu" stupids_julia!(a,b,e)







