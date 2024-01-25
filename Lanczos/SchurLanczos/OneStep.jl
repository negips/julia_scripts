#!/bin/julia


println("One Step Forward")

k        += 1

v         = A*v
w         = A'*w

Qlv,Rlv   = SchurLanczos.GetQRviews(Ql,Rl,k)
Qrv,Rrv   = SchurLanczos.GetQRviews(Qr,Rr,k)
γlv       = view(γl,1:k)
γrv       = view(γr,1:k)

SchurLanczos.BiOrthogonalize!(γlv,γrv,w,v,Qlv,Qrv,Rlv,Rrv,ngs)

hv        = view(H,:,k)
copyto!(hv,1,γlv,1,k)

bio = abs(w'*v)

if bio<BiOtol
  println("Vanishing BiOrthogonality: $(bio)")
else  
  αl,αr = SchurLanczos.BiOrthoScale!(w,v,x)
  if k<Lk
    SchurLanczos.UpdateQR!(Ql,Rl,w,k,ngs)
    SchurLanczos.UpdateQR!(Qr,Rr,v,k,ngs)
    hv[k+1] = αr
  end
  
end  

if (k==Lk)
  println("Subspace full: $Lk")
end  

println("Moved One Step: $k/$Lk")


  
