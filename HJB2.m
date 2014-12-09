function [Xi,v]=HJB2(t0,T,N,M1,M2,f0,F,g,U,Omega0,tol1,tol2)

[Omega,Lo,Up,MPsi,f0MPsi,FMPsi]=createsets2(Omega0,f0,F,N,M1,M2,U,t0,T,tol1)	
h=(T-t0)/N;
t=t0:h:T;
P=[1:-1/M1:0;0:1/M1:1];
Xi=Omega*P	;
v=zeros(N+1,M1+1);
v(N+1,:)=g(Xi(N+1,:));
Xi2=Xi(N+1,1);v2=v(N+1,1);
for j=1:M1
    [XiTmp,vXiTmp]=refine(Xi(N+1,j),Xi(N+1,j+1),v(N+1,j),v(N+1,j+1),g,tol2);	
    Xi2=[Xi2,XiTmp(2:end)];
    v2=[v2,vXiTmp(2:end)];
end
k=length(Xi2);
Xi(N+1,1:k+1)=[Xi2,nan];
v(N+1,1:k+1)=[v2,nan];
for i=N:-1:1
     #if(i<N)
     ind=isnan(Xi(i+1,:)); k=find(ind==1)-1;
     #else k=M1+1; end
     Xiipo=Xi(i+1,1:k);
     vXiipo=v(i+1,1:k);
     for j=1:M1+1
         v(i,j)=optimization(Xiipo,vXiipo,[Lo(i,j),Up(i,j)],i,j);
     end
     ## adaptive refinement of Xi
     ind=isnan(MPsi(i,:)); k=find(ind==1)-1;
     Psi=MPsi(i,1:k); f0Psi=f0MPsi(i,1:k); FPsi=FMPsi(i,1:k);
     Xi2=Xi(i,1);
     v2=v(i,1);
     for j=1:M1
         [XiTmp,vXiTmp]=refineV(Xi(i,j),Xi(i,j+1),v(i,j),v(i,j+1),Psi,f0Psi,FPsi,Xiipo,vXiipo,h,tol2,@(x) f0(t(i),x),@(x) F(t(i),x),U);
         Xi2=[Xi2,XiTmp(2:end)];
         v2=[v2,vXiTmp(2:end)];
     end
     k=length(Xi2);
     Xi(i,1:k+1)=[Xi2,nan];
     v(i,1:k+1)=[v2,nan];
end
v2=v;
for i=1:N+1
  k=find(isnan(Xi(i,:))); v(i,k:end)=nan; Xi(i,k:end)=nan;	
end
Xi=Xi(:,1:end-1);
v=v(:,1:end-1);
#v2(find(v2==0))=nan;
t=ones(size(Xi)(2),1)*(t0:(T-t0)/N:T);
surf(t',Xi,v)
			
end