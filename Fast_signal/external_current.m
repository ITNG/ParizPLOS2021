function Iext=external_current(N,range,dh,r,taud)
xrand=rand(N,range);
for ii=1:N
    X(ii,:)=xrand(ii,:)<r*dh/1000;
end

S=zeros(N,1);
Iext=zeros(N,range);
for ii=1:range
    fired=X(:,ii)==1;
    S(fired)=S(fired)+1;
    S=S-dh*S/taud;
    Iext(:,ii)=S;
end
    