function [Omega,Lo,Up,MPsi,f0MPsi,FMPsi]=createsets2(Omega0,f0,F,N,M1,M2,U,t0,T,tol)

## N...number of time steps
## M1..number of space steps for Evaluation of v (later refined)
## M2..number of space steps for Evaluation of f (later refined)

## tol..tolerance parameter for the adaptive discretization of f. 

Omega=zeros(N,2);
Lo=zeros(N,M1+1);
Up=zeros(N,M1+1);
h=(T-t0)/N;
t=t0:h:T;
Omega(1,:)=Omega0;
MPsi=zeros(N,M2+1);
f0MPsi=MPsi;
FMPsi=MPsi;
for i=1:N
    # calculate Omega_{i+1/2}
    k=(Omega(i,2)-Omega(i,1))/M1;
    x=Omega(i,1):k:Omega(i,2);
    for j=1:M1+1
       I=sort(x(j)+h*(f0(t(i),x(j))+F(t(i),x(j))*U));
       Lo(i,j)=I(1);
       Up(i,j)=I(2);
    end
    Psi0=min(Lo(i,:))-.1;
    Psi1=max(Up(i,:))+.1;
    # Ealuate f on an equidistant mesh
    k2=(Psi1-Psi0)/M2;
    Psi=Psi0:k2:Psi1;
    FPsi=F(t(i+1),Psi);
    f0Psi=f0(t(i+1),Psi);
    
    ### adaptive refinement for f
    Psi2=Psi(1); f0Psi2=f0Psi(1); FPsi2=FPsi(1);
    for j=1:M2
        [Psitmp,f0Psitmp,FPsitmp]=refine2(Psi(j),Psi(j+1),f0Psi(j),f0Psi(j+1),FPsi(j),FPsi(j+1),@(x) f0(t(i+1),x),@(x) F(t(i+1),x),tol,U);
        Psi2=[Psi2,Psitmp(2:end)]; f0Psi2=[f0Psi2,f0Psitmp(2:end)]; FPsi2=[FPsi2,FPsitmp(2:end)];
    end
    ## store the results of adaptive refinement of f for the adaptive refinement of v
    Psi=Psi2; f0Psi=f0Psi2; FPsi=FPsi2;
    M2tmp=length(Psi);
    MPsi(i,1:M2tmp+1)=[Psi,nan];
    f0MPsi(i,1:M2tmp+1)=[f0Psi,nan];
    FMPsi(i,1:M2tmp+1)=[FPsi,nan];
    ## calculate the reachable sets on an equidistant mesh
    for j=1:M1+1
        I=reachableset(x(j),U,h,Psi,f0Psi,FPsi,f0(t(i),x(j)),F(t(i),x(j)));
        Lo(i,j)=I(1);
        Up(i,j)=I(2);
    end
    Omega(i+1,1)=min(Lo(i,:));
    Omega(i+1,2)=max(Up(i,:));
end
figure(1)
plot(t,Omega(:,1),t,Omega(:,2));