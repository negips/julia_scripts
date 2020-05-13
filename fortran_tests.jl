println("Calling fortran subroutines")

using LinearAlgebra

n=5;
z=zeros(Float64,n)
w=zeros(Float64,n)

ccall((:zwgl_, "./speclib.so"), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Int32}), z,w,n)
println(z)

x = Float64[2.0];
y = Float64[0.];

#xp = Ref{x};
#yp = Ref{y};

ccall((:test_, "./fortran_sub.so"), Cvoid, (Ref{Float64}, Ref{Float64}), x, y)
println(y[1])

y2 = ccall((:test_, "./fortran_fcn.so"), Float64, (Ref{Float64},), x)
println(y2)

# DX = range(0.,1.0,length=10);
# DY = DX;
# n = length(DX)
# incx = incy = 1
# product = ccall((:ddot_, "libLAPACK"), Float64, (Ref{Int32}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Int32}), n, DX, incx, DY, incy)
# 
