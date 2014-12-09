function [Xi,vXi]=refineV(Xi1,Xi2,v1,v2,Psi,f0Psi,FPsi,Xiipo,vXiipo,h,tol2,f0,F,U)
XiNew=(Xi1+Xi2)/2;
f0New=f0(XiNew); FNew=F(XiNew);
INew=reachableset(XiNew,U,h,Psi,f0Psi,FPsi,f0New,FNew);
vNew=optimization(Xiipo,vXiipo,INew);
R=vNew-(v1+v2)/2;
if abs(R)>tol2	
   [Xil,vl]=refineV(Xi1,XiNew,v1,vNew,Psi,f0Psi,FPsi,Xiipo,vXiipo,h,tol2,f0,F,U);
   [Xir,vr]=refineV(XiNew,Xi2,vNew,v2,Psi,f0Psi,FPsi,Xiipo,vXiipo,h,tol2,f0,F,U);
   Xi=[Xil(1:end-1),Xir];
   vXi=[vl(1:end-1),vr];
else
    Xi=[Xi1,XiNew,Xi2]; vXi=[v1,vNew,v2];
end

end