function [Psi,fPsi]=refine(Psi1,Psi2,f1,f2,f,tol)
PsiNew=(Psi1+Psi2)/2; fNew=f(PsiNew);
R=fNew-(f1+f2)/2;
if abs(R)>tol	
   [Psil,fl]=refine(Psi1,PsiNew,f1,fNew,f,tol);
   [Psir,fr]=refine(PsiNew,Psi2,fNew,f2,f,tol);
   Psi=[Psil(1:end-1),Psir];
   fPsi=[fl(1:end-1),fr];
else
    Psi=[Psi1,PsiNew,Psi2]; fPsi=[f1,fNew,f2];
end
#Psi=Psi(2:end);
#fPsi=fPsi(2:end)
end