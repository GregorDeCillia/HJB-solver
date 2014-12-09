function [Xi,v]=HJB(t0,T,N,M1,M2,f0,F,g,U,Omega0)

[Omega,Lo,Up]=createsets(Omega0,f0,F,N,M1,M2,U,t0,T);	
P=[1:-1/M1:0;0:1/M1:1];
Xi=Omega*P;	
v=zeros(N+1,M1+1);
v(N+1,:)=g(Xi(N+1,:));
for i=N:-1:1
     vXi=v(i+1,:);
     for j=1:M1+1
     	i;j;
         v(i,j)=optimization(Xi(i+1,:),vXi,[Lo(i,j),Up(i,j)],i,j);
     end
end
t=ones(M1+1,1)*(t0:(T-t0)/N:T);
surf(t',Xi,v)
			
end