function [Us,Phi2] = Us(x,f0,F,U,h,Psi)
I2=x+h*(f0+F*U);
I2=sort(I2);
ind = (I2(1)<=Psi) .* (Psi<=I2(2));
Us=(Psi(find(ind==1))-x-h*f0)/(h*F);
find1=find(ind==1);
Phi2=[min(find1),max(find1)];
end 