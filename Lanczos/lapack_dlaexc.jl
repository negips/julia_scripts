function dlaexc!(wantq::BlasInt, j1::BlasInt, n1::BlasInt, n2::BlasInt, T::AbstractMatrix{Float64}, Q::AbstractMatrix{Float64})

            #BlasInt = LinearAlgebra.BlasInt

            #BlasInt = Int64

            LinearAlgebra.chkstride1(T, Q)
            n = LinearAlgebra.checksquare(T)
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            work = Vector{Float64}(undef, n)
            info = Ref{BlasInt}()
            #const libsrc = "/usr/lib/julia/libblastrampoline.so.5"
            ccall((:dlaexc_, localsrc), Cvoid,
                  (Ref{BlasInt},Ref{BlasInt},
                   Ptr{Float64},Ref{BlasInt},Ptr{Float64},Ref{BlasInt},
                   Ref{BlasInt},Ref{BlasInt},Ref{BlasInt},
                   Ptr{Float64},Ptr{BlasInt},Clong),
                   wantq,n,
                   T,ldt,Q,ldq,
                   j1,n1,n2,
                   work,info,1)


                 # (Ref{UInt8},  Ref{BlasInt},
                 #  Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 #  Ref{BlasInt}, Ref{BlasInt},
                 #  Ptr{$elty}, Ptr{BlasInt}, Clong),
                 # compq, n,
                 # T, ldt, Q, ldq,
                 # ifst, ilst,
                 # work, info, 1)
            #chklapackerror(info[])
            T, Q
        end




#for (trexc, trsen, tgsen, elty) in
#    ((:dtrexc_, :dtrsen_, :dtgsen_, :Float64),
#     (:strexc_, :strsen_, :stgsen_, :Float32))
#    @eval begin
#        # *     .. Scalar Arguments ..
#        #       CHARACTER          COMPQ
#        #       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
#        # *     ..
#        # *     .. Array Arguments ..
#        #       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
#
#        function trexc!(compq::AbstractChar, ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty})
#            chkstride1(T, Q)
#            n = checksquare(T)
#            ldt = max(1, stride(T, 2))
#            ldq = max(1, stride(Q, 2))
#            work = Vector{$elty}(undef, n)
#            info = Ref{BlasInt}()
#            ccall((@blasfunc($trexc), libblastrampoline), Cvoid,
#                  (Ref{UInt8},  Ref{BlasInt},
#                   Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
#                   Ref{BlasInt}, Ref{BlasInt},
#                   Ptr{$elty}, Ptr{BlasInt}, Clong),
#                  compq, n,
#                  T, ldt, Q, ldq,
#                  ifst, ilst,
#                  work, info, 1)
#            chklapackerror(info[])
#            T, Q
#        end

