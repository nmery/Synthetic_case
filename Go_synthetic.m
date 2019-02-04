seed=821234651;
rng(seed);
load tab.mat;
modelg=[1 1;4 20];cg=[0.3;0.7];
yz=tab(:,:,1);
[model,c]=transcov(modelg,cg,yz);
n=10000;
vn=[3 6 9 12 24 36];  % number of samples
vnp=[1 2 3 4 5 6];    % id boreholes
lse=20;     % external window size
lsi=0;      % internal window size
lc=3;       % length (composite)

% cas=1 simulated 5 boreholes per quadrant
% cas=2 simulated 3 boreholes per quadrant

for cas=1:2
    
    for ii=1:length(vn)
        m=vn(ii);
        mp=vnp(ii);
        
        if cas==1
            [be,bt,osrt,osre,ke,be_cr,bt_cr,osrt_cr,osre_cr,ke_cr,stat,zv,ze,ze_cr]=test_borehole_5(n,modelg,cg,model,c,lse,m,mp,yz,lc,lsi,cas);
            statall{ii}=stat;
        elseif cas==2
            [be,bt,osrt,osre,ke,be_cr,bt_cr,osrt_cr,osre_cr,ke_cr,stat,zv,ze,ze_cr]=test_borehole_3(n,modelg,cg,model,c,lse,m,mp,yz,lc,lsi,cas);
            statall{ii}=stat;
        end
        
        vbe(ii)=be;
        vbt(ii)=bt;
        vosrt(ii)=osrt;
        vosre(ii)=osre;
        vke(ii)=ke;
        
        vbe_cr(ii)=be_cr;
        vbt_cr(ii)=bt_cr;
        vosrt_cr(ii)=osrt_cr;
        vosre_cr(ii)=osre_cr;
        vke_cr(ii)=ke_cr;
    end
    
    
    figure(10+cas);clf
    plot(vnp(3:end),vbt(3:end),'-b','linewidth',2);
    hold on
    plot(vnp(3:end),vbe(3:end),'--b','linewidth',2);
    plot(vnp(3:end),vosrt(3:end),'-r',vnp(3:end),vosre(3:end),'--r','linewidth',2);
    plot(vnp(3:end),vke(3:end),'-k','linewidth',2)
    axis([2.7 6.5 -7 100])
    grid on
    hold off
    h=legend('$SR th$','$SR exp$','$OSR th$','$OSR exp$','$KE$');
    set(h,'Interpreter','Latex','Location','southeast');
    set(gca,'TickLabelInterpreter','latex')
    box off
    xticks([3 4 5 6])
    xticklabels({'3','1xQ','2xQ','3xQ'})
    xlabel('Number of boreholes ','Interpreter','Latex','FontSize', 12)
    ylabel('$\%$','Interpreter','Latex','FontSize', 12)
    print(10+cas,'-depsc2',['KNA_OK_cas' num2str(cas) '.eps']);
    close(figure(10+cas));
    
    figure(11+cas);clf;
    plot(vnp(3:end),vbt_cr(3:end),'-b','linewidth',2);
    hold on
    plot(vnp(3:end),vbe_cr(3:end),'--b','linewidth',2);
    plot(vnp(3:end),vosrt_cr(3:end),'-r',vnp(3:end),vosre_cr(3:end),'--r','linewidth',2);
    plot(vnp(3:end),vke_cr(3:end),'-k','linewidth',2)
    h=legend('$SR th$','$SR exp$','$OSR th$','$OSR exp$','$KE$');
    set(h,'Interpreter','Latex','Location','southeast');
    set(gca,'TickLabelInterpreter','latex')
    box off
    xticks([3 4 5 6])
    xticklabels({'3','1xQ','2xQ','3xQ'})
    xlabel('Number of boreholes','Interpreter','Latex','FontSize', 12)
    ylabel('$\%$','Interpreter','Latex','FontSize', 12)
    axis([2.7 6.5 -7 100])
    grid on
    hold off
    print(11+cas,'-depsc2',['KNA_CK_cas' num2str(cas) '.eps']);
    close(figure(11+cas));
    
end