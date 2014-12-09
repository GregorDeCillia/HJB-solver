function I = reachableset(x,U,h,Psi,f0Psi,FPsi,f0,F)
# Calculate the subintervals of U with the correpondig Psi break points
if (F==0)
  I=x+h*f0;
  ind=Psi<=I;
  j=max(find(ind==1));
  lambda=(I-Psi(j))/(Psi(j+1)-Psi(j));
  fPsiU=(1-lambda).*(f0Psi(j)+FPsi(j)*U)+lambda.*(f0Psi(j+1)+FPsi(j+1)*U);
  Phih=f0+F*U+fPsiU;
  Psi=[0,inf];
  n=0;

else
	
[Us,Phi2] = Us(x,f0,F,U,h,Psi);
# Initialize for the loop
if isempty(Us)     #     ???????
  I2=x+h*(f0+F*U);
  I=sort(I2);	
  ind=I2(1)>=Psi;
  j=max(find(ind==1));
  Umm=(Psi(j)-x-h*f0)/(h*F);
  Upp=(Psi(j+1)-x-h*f0)/(h*F);
  lambda=(U-Umm)/(Upp-Umm);
  fPsiU=(1-lambda).*(f0Psi(j)+FPsi(j)*U)+lambda.*(f0Psi(j+1)+FPsi(j+1)*U);
  Phih=f0+F*U+fPsiU;
  Psi=(1-lambda)*Psi(j)+lambda*Psi(j+1);
  n=0;
else

  n=length(Us);
  j=Phi2(1)-1;
  fPsiU=zeros(1,n);
  Phih=fPsiU;
  for i=1:n
     fPsiU(i)=f0Psi(j+i)+FPsi(j+i)*Us(i);
     Phih(i)=f0+F*Us(i)+fPsiU(i);
  end
  ### wiese jedem randwert von Us einen Randwert von U zu
  if (F>0), U1=U(1); U2=U(2);
  else U1=U(2); U2=U(1);
  end
  ### Berechne PsiU, Phih und Psi für die Randwerte
  Umm=(Psi(j)-x-h*f0)/(h*F);
  lambda=(U1-Umm)/(Us(1)-Umm);
  fPsiUm=(1-lambda)*(f0Psi(j)+FPsi(j)*U1)+lambda*(f0Psi(j+1)+FPsi(j+1)*U1);
  Phihm=f0+F*U1+fPsiUm;
  Psim=(1-lambda)*Psi(j)+lambda*Psi(j+1);

  Upp=(Psi(j+n+1)-x-h*f0)/(h*F);
  lambda=(U2-Upp)/(Us(n)-Upp);
  fPsiUp=(1-  lambda)*(f0Psi(j+n+1)+FPsi(j+n+1)*U2)+lambda*(f0Psi(j+n)+FPsi(j+n)*U2);
  Phihp=f0+F*U2+fPsiUp;
  Psip=(1-lambda)*Psi(j+n+1)+lambda*Psi(j+n);

  fPsiU=[fPsiUm,fPsiU,fPsiUp];
  Phih=[Phihm,Phih,Phihp];
  Psi=[Psim,Psi(Phi2(1):Phi2(2)),Psip];
end

end
inn=Phih;
for i=1:n+1
     dlambda=h*F/(Psi(i+1)-Psi(i));
     # calculate the derivatives of (u->f+fPsi) on the boundaries of the subintervals of U
     d1=dlambda*(fPsiU(i+1)-fPsiU(i))+F;
     d2=d1+FPsi(i+1);
     d1=d1+FPsi(i);
     if (d1*d2<0) # necessary and suficiient condition for inner optimum
          lambda=-d1/(d2-d1);
          inn(i)=Phih(i)+(Psi(i+1)-Psi(i))*(d1*lambda+(d2-d1)*lambda^2/2); # taylor formula 
     end
end
inn;
Phih;

inn==Phih;

I=[min([inn,Phih]),max([inn,Phih])];
I=x+h/2*I;
end