function [be,bt,osrt,osre,ke,be_cr,bt_cr,osrt_cr,osre_cr,ke_cr,stat,zv,ze,ze_cr]=test_borehole_3(n,modelg,cg,model,c,lse,m,mp,yz,lc,lsi,cas)

x0=grille3(-2,2,1,-2,2,1,-2,2,1);
k0=covardm(x0,x0,model,c);
sv=mean(k0(:));
nn=size(x0,1)+m;
zv=zeros(n,1);
ze=zeros(n,1);
ze_cr=ze;
bt1=zeros(n,1);
bt1_cr=bt1;
osrt=zeros(n,1);
osrt_cr=osrt;
s=zeros(n,1);
KE=zeros(n,1);
KE_cr=KE;
index=zeros(n,1);
nbimag=0;

for i=1:n
    bh1x=(-lsi+lse).*rand(1,1)-lse;
    bh1y=(lse-lsi).*rand(1,1)+lsi;
    bh1x2=(-lsi+lse).*rand(1,1)-lse;
    bh1y2=(lse-lsi).*rand(1,1)+lsi;
    bh1x3=(-lsi+lse).*rand(1,1)-lse;
    bh1y3=(lse-lsi).*rand(1,1)+lsi;
    bh1z=-(0.5+rand)*lc:lc:1.5*lc;
    bh2x=(lse-lsi).*rand(1,1)+lsi;
    bh2y=(lse-lsi).*rand(1,1)+lsi;
    bh2x2=(lse-lsi).*rand(1,1)+lsi;
    bh2y2=(lse-lsi).*rand(1,1)+lsi;
    bh2x3=(lse-lsi).*rand(1,1)+lsi;
    bh2y3=(lse-lsi).*rand(1,1)+lsi;
    bh2z=-(0.5+rand)*lc:lc:1.5*lc;
    bh3x=(-lsi+lse).*rand(1,1)-lse;
    bh3y=(-lsi+lse).*rand(1,1)-lse;
    bh3x3=(-lsi+lse).*rand(1,1)-lse;
    bh3y3=(-lsi+lse).*rand(1,1)-lse;
    bh3x2=(-lsi+lse).*rand(1,1)-lse;
    bh3y2=(-lsi+lse).*rand(1,1)-lse;
    bh3z=-(0.5+rand)*lc:lc:1.5*lc;
    bh4x=(lse-lsi).*rand(1,1)+lsi;
    bh4y=(-lsi+lse).*rand(1,1)-lse;
    bh4x2=(lse-lsi).*rand(1,1)+lsi;
    bh4y2=(-lsi+lse).*rand(1,1)-lse;
    bh4x3=(lse-lsi).*rand(1,1)+lsi;
    bh4y3=(-lsi+lse).*rand(1,1)-lse;
    bh4z=-(0.5+rand)*lc:lc:1.5*lc;
    bh1=[[bh1x;bh1x;bh1x],[bh1y;bh1y;bh1y],bh1z'];
    bh2=[[bh2x;bh2x;bh2x],[bh2y;bh2y;bh2y],bh2z'];
    bh3=[[bh3x;bh3x;bh3x],[bh3y;bh3y;bh3y],bh3z'];
    bh4=[[bh4x;bh4x;bh4x],[bh4y;bh4y;bh4y],bh4z'];
    bh5=[[bh1x2;bh1x2;bh1x2],[bh1y2;bh1y2;bh1y2],bh1z'];
    bh6=[[bh2x2;bh2x2;bh2x2],[bh2y2;bh2y2;bh2y2],bh2z'];
    bh7=[[bh3x2;bh3x2;bh3x2],[bh3y2;bh3y2;bh3y2],bh3z'];
    bh8=[[bh4x2;bh4x2;bh4x2],[bh4y2;bh4y2;bh4y2],bh4z'];
    bh9=[[bh1x3;bh1x3;bh1x3],[bh1y3;bh1y3;bh1y3],bh1z'];
    bh10=[[bh2x3;bh2x3;bh2x3],[bh2y3;bh2y3;bh2y3],bh2z'];
    bh11=[[bh3x3;bh3x3;bh3x3],[bh3y3;bh3y3;bh3y3],bh3z'];
    bh12=[[bh4x3;bh4x3;bh4x3],[bh4y3;bh4y3;bh4y3],bh4z'];
    x2=[bh1;bh5;bh9;bh2;bh6;bh10;bh3;bh7;bh11;bh4;bh8;bh12];
    
    x1=[0,0];
    [n1,~]=size(x1);
    [n2,~]=size(x2);
    
    h=sqrt((x1(:,1)*ones(1,n2)-ones(n1,1)*x2(:,1)').^2+(x1(:,2)*ones(1,n2)-ones(n1,1)*x2(:,2)').^2);
    h=h';
    x21=x2(1:9,:);
    x22=x2(10:18,:);
    x23=x2(19:27,:);
    x24=x2(28:36,:);
    
    [~,id0]=sort(h(:,:));
    [~,id1]=sort(h(1:9,:));
    [~,id2]=sort(h(10:18,:));
    [~,id3]=sort(h(19:27,:));
    [~,id4]=sort(h(28:36,:));
    
    if mp==1
        x=x2(id0(1:3),:);
    elseif mp==2
        x=x2(id0(1:6),:);
    elseif mp==3
        x=x2(id0(1:9),:);
    elseif mp==4
        x=[x21(id1(1:3),:);x22(id2(1:3),:);x23(id3(1:3),:);x24(id4(1:3),:)];
    elseif mp==5
        x=[x21(id1(1:6),:);x22(id2(1:6),:);x23(id3(1:6),:);x24(id4(1:6),:)];
    elseif mp==6
        x=x2;
    end
    
    kx=covardm(x,x,model,c);
    kxx0=covardm(x,x0,model,c);
    k=[kx kxx0;kxx0' k0];
    kg=[kx ones(m,1);ones(1,m) 0];
    kd=[mean(kxx0,2);1];
    
    l=kg\kd;
    ll=l(1:m);
    varze=ll'*kx*ll;
    t=diag(l'*kd);
    s(i)=sv-t';
    
    %KNA
    bt1(i)=100*ll'*kd(1:m)/varze;
    osrt(i)=100*((sv-varze)/sv);
    KE(i)=100.*(sv-s(i))./sv;
    
    % Cressie
    ki=inv(kx);
    ski=sum(ki(:));
    l=ki*kd(1:m);
    k0l=kd(1:m)'*l;
    mu2=sqrt(l'*kd(1:m)*ski-sum(l)^2)/sqrt(ski*sv-1);
    mu1=-mu2/ski+sum(l)/ski;
    lc=1/mu2*kd(1:m)'*ki-mu1/mu2*sum(ki);
    lc=lc';
    if sum(abs(imag(lc)))>0
        s_cr(i)=NaN;
        bt1_cr(i)=NaN;
        osrt_cr(i)=NaN;
        KE_cr(i)=NaN;
        index(i)=0;
        varze_cr=NaN;
        lc=ll;
        nbimag=nbimag+1;
    else
        varze_cr=lc'*kx*lc;
        s_cr(i)=sv+varze_cr-2*lc'*kd(1:m);
        bt1_cr(i)=100*lc'*kd(1:m)/varze_cr;
        osrt_cr(i)=100*((sv-varze_cr)/sv);
        KE_cr(i)=100.*(sv-s_cr(i))./sv;
        index(i)=i;
    end
    
    % simulation of block and sample values
    kx=covardm(x,x,modelg,cg);
    kxx0=covardm(x,x0,modelg,cg);
    k0=covardm(x0,x0,modelg,cg);
    k=[kx kxx0;kxx0' k0];
    
    u=chol(k);
    y=u'*randn(nn,1);
    z=interp1(yz(:,1),yz(:,2),y,'linear','extrap');
    zv(i)=mean(z(m+1:end));
    ze_cr(i)=lc'*z(1:m);
    ze(i)=ll'*z(1:m);
    
end

[(sum(ze<0))/n (sum(ze_cr)<0)/n nbimag/n];

%remove negative values
ze_cr(ze_cr<0)=0;
ze(ze<0)=0;
index(index==0)=[];

[mean(zv) mean(ze) mean(ze_cr)];

be=regres(zv,ze,1);
be=be(1)*100;
bt=nanmean(bt1);
osrt=nanmean(osrt);
ke=nanmean(KE);
osre=100*((var(zv)-var(ze))/var(zv));

% cressie
be_cr=regres(zv(index),ze_cr(index),1);
be_cr=be_cr(1)*100;
if be_cr==0
    be_cr=NaN;
else
    be_cr=be_cr;
end
bt_cr=nanmean(bt1_cr);
osrt_cr=nanmean(osrt_cr);
ke_cr=nanmean(KE_cr);
osre_cr=100*((var(zv(index))-var(ze_cr(index)))/var(zv(index)));

[osrt osre osrt_cr osre_cr];

% compute recovery functions
p=[0.025:0.05:1]';
c=quantile(zv,p);
for i=1:length(c)
    tc(i)=mean(zv>=c(i));
    tce(i)=mean(ze>=c(i));
    tc_cr(i)=mean(ze_cr>=c(i));
    
    mc(i)=mean(zv(zv>=c(i)));
    mce(i)=mean(ze(ze>=c(i)));
    mc_cr(i)=mean(ze_cr(ze_cr>=c(i)));
    
    pc(i)=tc(i)*(mc(i)-c(i));
    pce(i)=tce(i)*(mce(i)-c(i));
    pc_cr(i)=tc_cr(i)*(mc_cr(i)-c(i));
    
end
stat.propneg=nanmean(ze<0);
stat.propneg_cr=nanmean(ze_cr<0);
stat.prop_imag=nbimag/n;
stat.n=n;
stat.nbneg=sum(ze<0);
stat.nbneg_cr=sum(ze_cr<0);
stat.nbimag=nbimag;
stat.rmse=sqrt(nanmean((zv-ze).^2));
stat.rmse_cr=sqrt(nanmean((zv-ze_cr).^2));
t=corrcoef(zv,ze);
stat.corr=t(1,2);
tcc=corrcoef(zv,ze_cr);
stat.corr_cr=tcc(1,2);
stat.tc=tc;
stat.mc=mc;
stat.pc=pc;
stat.tce=tce;
stat.mce=mce;
stat.pce=pce;
stat.tc_cr=tc_cr;
stat.mc_cr=mc_cr;
stat.pc_cr=pc_cr;

[sqrt(nanmean((zv-ze).^2)) sqrt(nanmean((zv-ze_cr).^2)) t(1,2) tcc(1,2)];


figure(100+m);clf
plot(c,tc,'--r',c,tce,'-b',c,tc_cr,'-k','linewidth',2)
title(['Neighbors ',num2str(m)],'fontsize',15)
grid on
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','northeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Tonnage','Interpreter','Latex','FontSize', 15)
title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)
print(100+m,'-depsc2',['CCT' num2str(mp) '_cas' num2str(cas) '.eps']);
close(figure(100+m));


figure(101+m)
plot(c,mc,'--r',c,mce,'-b',c,mc_cr,'-k','linewidth',2)
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','southeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title(['Curve ore grade - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)
grid on
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Ore grade','Interpreter','Latex','FontSize', 15)
print(101+m,'-depsc2',['CCG' num2str(mp) '_cas' num2str(cas) '.eps']);
close(figure(101+m));

figure(102+m)
plot(c,pc,'--r',c,pce,'-b',c,pc_cr,'-k','linewidth',2)
h=legend('Real','OK','CK');
set(h,'Interpreter','Latex','Location','northeast','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title(['Curve Conv. Profit - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 15)
grid on
box off
xlabel('Cut-off','Interpreter','Latex','FontSize', 15)
ylabel('Conventional profit','Interpreter','Latex','FontSize', 15)
print(102+m,'-depsc2',['CCC' num2str(mp) '_cas' num2str(cas) '.eps']);
close(figure(102+m));

drawnow
