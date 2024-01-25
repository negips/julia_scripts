#!/bin/julia


println("Restart Step: Nev = $Nev")

Sr = schur(H)

select = SchurLanczos.SelectEigenvalues(Sr.values,Nev)

T = copy(Sr.Schur)
Z = Matrix(diagm(ones(vt,Lk)))

LinearAlgebra.LAPACK.trsen!(select,T,Z)
# Z is the transofrmation such that the new Schur factorization Tnew = Z'*T*Z

Znew = Z*Sr.Z

RrZ  = Rr*Znew[:,1:Nev]

qrz  = qr(RrZ)

println("Done.")

#SchurLanczos.BiOrthogonalize!(γlv,γrv,w,v,Qlv,Qrv,Rlv,Rrv,ngs)
#
#hv        = view(H,:,k)
#copyto!(hv,1,γlv,1,k)
#
#bio = abs(w'*v)
#
#if bio<BiOtol
#  println("Vanishing BiOrthogonality: $(bio)")
#else  
#  αl,αr = SchurLanczos.BiOrthoScale!(w,v,x)
#  if k<Lk
#    SchurLanczos.UpdateQR!(Ql,Rl,w,k,ngs)
#    SchurLanczos.UpdateQR!(Qr,Rr,v,k,ngs)
#    hv[k+1] = αr
#  end
#  
#end  
#
#if (k==Lk)
#  println("Subspace full: $Lk")
#end  
#
#println("Moved One Step: $k/$Lk")


  
