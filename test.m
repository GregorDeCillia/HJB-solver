f0=@(t,x) sin(x)
F=@(t,x) cos(x)
g=@(x) sin(2*x)
Omega0=[0,1]
t0=1, T=2
N=10
M1=20
saver=zeros(10,M1+1);
for (i=1:10)
     M2=2^i
     [Xi,v]=HJB(t0,T,N,M1,M2,f0,F,g,U,Omega0);
     saver(i,:)=v(1,:)
end
err=abs(saver-ones(10,1)*saver(end,:))

loglog(2.^(1:9),err(1:end-1,:))
### shows the error is quadratic w.r.t. M2