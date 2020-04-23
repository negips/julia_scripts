using ImageView, Images, ImageMagick, Printf, FileIO, LinearAlgebra# , Arpack

basefilename="market_";
m = 2; # Numer of frames
fname = @sprintf("%s%04d.png",basefilename,1);

println(fname)
img = load(fname);
imshow(img)

sz = size(img); szv=sz[1]*sz[2];

R = float(red.(img));
G = float(green.(img));
B = float(blue.(img));
v = [reshape(R,szv,1);reshape(G,szv,1);reshape(B,szv,1)]

A = zeros(size(v,1),m);
p,q = size(A);
for k=1:m
    fname = @sprintf("%s%04d.png",basefilename,k);
    println(fname)
    img = load(fname);
    
    sz = size(img); szv=sz[1]*sz[2];
    
    R = float(red.(img));
    G = float(green.(img));
    B = float(blue.(img));
    v = [reshape(R,szv,1);reshape(G,szv,1);reshape(B,szv,1)]
    A[:,k] = v;
     
end

B = zeros(p+q,p+q);

B = [zeros(p,p A'; A zeros(q,q)];

# Define B and C and use arnoldi to get its eigenvalues

#U,S,V = svd(A); S=diagm(0=>S);

#A_1=U[:,1]*S[1,1]*V[:,1]';

#v=A_1[:,1]; # Third column of A = third frame
#vv=reshape(v,szv,3);
#R=reshape(vv[:,1],sz[1],sz[2]);
#G=reshape(vv[:,2],sz[1],sz[2]);
#B=reshape(vv[:,3],sz[1],sz[2]);
#newimg=RGB.(R,G,B)
#imshow(newimg)
#fname = "india_driving_bestrank.png";
#save(fname,newimg)

# The reverse of the problem

#for k in 1:m
#    v=A[:,k]; # Third column of A = third frame
#    vv=reshape(v,szv,3);
#    R=reshape(vv[:,1],sz[1],sz[2]);
#    G=reshape(vv[:,2],sz[1],sz[2]);
#    B=reshape(vv[:,3],sz[1],sz[2]);
#    newimg=RGB.(R,G,B)
#    imshow(newimg)
#end
