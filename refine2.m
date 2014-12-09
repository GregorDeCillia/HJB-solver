function [Psi,f0Psi,FPsi]=refine2(Psi1,Psi2,f01,f02,F1,F2,f0,F,tol,U)
PsiNew=(Psi1+Psi2)/2; f0New=f0(PsiNew); FNew=F(PsiNew);
R1=f0New-(f01+f02)/2; R2=FNew-(F1+F2)/2;
if max(abs(R1+U.*R2))>tol	
   [Psil,f0l,Fl]=refine2(Psi1,PsiNew,f01,f0New,F1,FNew,f0,F,tol,U);
   [Psir,f0r,Fr]=refine2(PsiNew,Psi2,f0New,f02,FNew,F2,f0,F,tol,U);
   Psi=[Psil(1:end-1),Psir];
   f0Psi=[f0l(1:end-1),f0r];
   FPsi=[Fl(1:end-1),Fr];
else
    Psi=[Psi1,PsiNew,Psi2]; f0Psi=[f01,f0New,f02]; FPsi=[F1,FNew,F2];
end
#Psi=Psi(2:end);
#fPsi=fPsi(2:end)
end