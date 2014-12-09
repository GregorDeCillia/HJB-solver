function v = optimization(Xi,vXi,I,i,j)
	i;j;
ind=(Xi<=I(2)).*(I(1)<=Xi);
if (sum(ind)==0)
	i,j,I,Xi,vXi,I
   ind=(Xi<=I(2));
   j=max(find(ind==1));
   Xi(j)
   lambda=(I-Xi(j))/(Xi(j+1)-Xi(j))
   v=(1-lambda).*vXi(j)+lambda*vXi(j+1);
   v=min(v) 	
else
find1=find(ind==1);
v=min(vXi(find1));
v=vXi(find1);
n=min(find1);
if n>1
lambda0=(I(1)-Xi(n-1))/(Xi(n)-Xi(n-1));
v0=(1-lambda0)*vXi(n-1)+lambda0*vXi(n);
v=[v,v0];
end
m=max(find1);
if m<length(Xi)
lambda1=(I(2)-Xi(m))/(Xi(m+1)-Xi(m));
v1=(1-lambda1)*vXi(m)+lambda1*vXi(m+1);
v=[v,v1];
end
#i,j,v
v=min(v);
end