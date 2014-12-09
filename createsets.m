function [Omega,Lo,Up]=createsets(Omega0,f0,F,N,M1,M2,U,t0,T)
## N...number of time steps
## M1..number of space steps for Evaluation of v
## M2..number of space steps for Evaluation of f

Omega=zeros(N,2);
Lo=zeros(N,M1+1);
Up=zeros(N,M1+1);
h=(T-t0)/N;
t=t0:h:T;
Omega(1,:)=Omega0;
for i=1:N
    k=(Omega(i,2)-Omega(i,1))/M1;
    x=Omega(i,1):k:Omega(i,2);
    for j=1:M1+1
       I=sort(x(j)+h*(f0(t(i),x(j))+F(t(i),x(j))*U));
       Lo(i,j)=I(1);
       Up(i,j)=I(2);
    end
    Psi0=min(Lo(i,:))-.1;
    Psi1=max(Up(i,:))+.1;
    k2=(Psi1-Psi0)/M2;
    Psi=Psi0:k2:Psi1;
    FPsi=F(t(i+1),Psi);
    f0Psi=f0(t(i+1),Psi);
    for j=1:M1+1
    	i,j,x(j),f0(t(i),x(j)),F(t(i),x(j)),Psi,f0Psi,FPsi
        I=reachableset(x(j),U,h,Psi,f0Psi,FPsi,f0(t(i),x(j)),F(t(i),x(j)));
        Lo(i,j)=I(1);
        Up(i,j)=I(2);
    end
    Omega(i+1,1)=min(Lo(i,:));
    Omega(i+1,2)=max(Up(i,:));
end
figure(1)
plot(t,Omega(:,1),t,Omega(:,2));