function [model,c]=transcov(modelg,cg,yz);
seed=915;
rng(seed);
yz=[-7 0;yz;7 1.5*max(yz(:,2))];

hmax=1.2*max(modelg(:,2));
dh=hmax/100;
h=[0, 0.00001:dh:hmax]';
k=covardm(0,h,modelg,cg);
n=100000;
y=randn(n,1);
z=interp1(yz(:,1),yz(:,2),y,'linear','extrap');
kt=zeros(size(k));
kt(1)=var(z);

for i=2:length(k);
    kk=[k(1) k(i);k(i) k(1)];
    u=chol(kk);
    yy=u'*randn(2,n);
    t1=interp1(yz(:,1),yz(:,2),yy(1,:)');
    t2=interp1(yz(:,1),yz(:,2),yy(2,:)');
    tt=cov(t1,t2);
    kt(i)=tt(1,2);
end

plot(h(:),kt(:),'-r','linewidth',2);
tt=kt(1)-kt(2:end);
gexp=[h(2:end) ones(length(h)-1,1) tt(:)];
[model,c,dif]=ajustermulti(gexp,[1 1;2 hmax/6;4 hmax/2],kt(1)*[1/3;1/3;1/3],[1 1;1 0;1 0],[0;0;0],[0 0],1,1,1)

