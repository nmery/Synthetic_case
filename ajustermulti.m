function [modela,ca,dif]=ajustermulti(gexp,model,c,modelcode,ccode,vdir,vreg,expoh,hmin);
%  fonction pour ajuster simultanément un modele dans plusieurs directions
%
%  syntaxe [modela,ca,dif]=ajustermulti(gexp,model,c,modelcode,ccode,vdir,vreg,expoh,hmin)
%
%  input   gexp: variog. experimental dans l'order: [h,n(h),g(h)]
%          model: modeles de variogramme (comme cokri) (attention ne pas mettre de portées 0 meme pour effet de pepite)
%          c:  vecteur des plateaux (comme cokri)
%          modelcode: ayant meme taille que model avec des 0 et des 1 (0 si libre, 1 si imposé), les parametres imposes
%                     ne sont pas optimisés. Le type, les portées et les rotations d'un effet de pépite ne sont pas modélisés
%          ccode: vecteur avec 0 ou 1 (0 libre et 1 imposé) pour l'ajustement des plateaux (c) des differentes structures
%          vdir  la matrice ndir x 2 direction et pendage ayant servis a calculer gexp (en 2D mettre pendage=0)
%          vreg: vecteur ndir x1 contenant les régularisations utilisées pour calculer gexp
%  output  modela: modele ajusté (comme model en input)       
%          ca: seuils ajustés;
%          dif: statistique d'ajustement
%
%  ATTENTION: Les angles de rotation de model sont donnés selon la convention cokri (trigonométrique)
%             en 3D, les rotations sont selon x, y et z dans model mais sont effectuées dans l'ordre z, y et x (antihoraire, syteme main-droite)
%  PAR CONTRE : Les angles de vdir sont des azimut et pendage (convention géologique; comme dans varioexp2d (azimut seulement) et varioexp3d (azimut et pendage)).
%
%  Auteur: D. Marcotte, juin 2002  

[n,p]=size(vdir);
if p==1
    vdir=[vdir zeros(n,1)];
end

u=poletocart(vdir);
ndir=size(vdir,1);
x=[];g=[];w=[];

for idir=1:ndir
    id=gexp(:,2,idir)>0;
    x=[x; gexp(id,1,idir)*u(idir,:)];
    g=[g; gexp(id,3,idir)];
    w=[w;gexp(id,2,idir).*ones(sum(id),1)./(hmin+gexp(id,1,idir)).^expoh];  % pondération pour nombre de paires et distance
end

% vérifier si l'on est avec un cas 2D anisotrope ou 3D anisotrope
[n,p]=size(model);
deuxd=0;

if p==4;
    id=model(:,2:3)==0;
    model(id)=0.0001;
    deuxd=1;
    model=[model(:,1:3),model(:,3),zeros(n,2),model(:,4)]; % compléter model pour rendre 3d
    modelcode=[modelcode(:,1:3),ones(n,3),modelcode(:,4)]; % les parametres ajoutés pour passer 3d sont fixés
elseif p==7
    id=model(:,2:4)==0;
    model(id)=0.0001;
    
end

id=model(:,1)==1;
model(id,2:end)=1;
modelcode(id,2:end)=1;
modelcode(:,1)=1;


options=optimset('TolX',0.00001,'TolFun',0.00001,'MaxFunEvals',13500);

id1=modelcode==0;p1=sum(sum(id1));
x0libre=reshape(model(id1),p1,1);
id2=ccode==0;p2=sum(sum(id2));
x0libre=[x0libre;c(id2)];

if ~isempty(x0libre)
    [t1,dif,flag]=fminsearch('gocovardm',x0libre,options,x,g,w,model,c,id1,id2);
    if flag==0
        display('Arret sur nombre d''itérations')
    else
        display('Arret a un minimum local')
    end
    
    if p1>0;
        t11=t1(1:p1);
        model(id1)=t11;
    end
    if p2>0
        t21=t1(p1+1:p1+p2);
        c(id2)=t21;
    end
else
   dif=gocovardm(x0libre,x,g,w,model,c,id1,id2);
   modela=model;
   ca=c;

end

if deuxd==1
    modela=model(:,[1:3 7]);
else
    modela=model;
end

ca=c;
visualiser(gexp,model,c,vdir,vreg,deuxd)


function visualiser(gexp,model,c,vdir,vreg,deuxd);
%
% syntaxe visualiser(gexp,model,c,vdir);;
%
% input   model et c comme cokri
%            gexp: variog. experimental de la variable (h, g(h)) 
%            vdir: azimuts et pendage a considérer
%            vreg: régularisations
%            deuxd: code pour indiquer 2d
%

% Définition de la grille

ndir=size(vdir,1);
nplot=ceil(sqrt(ndir));
if nplot*nplot==ndir
    nplot=nplot+1;
end

dh=[];hmax=[];gmax=[];


for idir=1:ndir
    id=gexp(:,2,idir)>0;
    dh=[dh;min(gexp(id,1,idir))];
    hmax=[hmax;max(gexp(id,1,idir))];
    gmax=[gmax;max(gexp(id,3,idir))];
end
dh=0.1*min(dh);
hmax=1.2*max(hmax);
gmax=1.2*max(gmax);

h=[dh:dh:hmax]';

% Boucler sur les directions
u=poletocart(vdir);

for idir=1:ndir
    nn=find(gexp(:,2,idir)>0);
    id=gexp(:,2,idir)>0;
    x=h*u(idir,:);
    k=covardm(x,[0,0,0],model,c);
    gt=sum(c)-k;
    %subplot(nplot,nplot,idir)
    set(gca,'TickLabelInterpreter','latex')
    box off
    plot(gexp(id,1,idir),gexp(id,3,idir),'+','markersize',5);hold on
    
    %for j=1:length(nn);
     %   t=text(gexp(nn(j),1,idir),gexp(nn(j),3,idir),[' ',num2str(gexp(nn(j),2,idir))],'fontsize',8);
      %  set(t,'fontsize',8,'tag','nombrepaires');
    %end
    
    




    plot(h,gt);hold off
    axis([0 hmax, 0, gmax]);
    set(gca,'TickLabelInterpreter','latex','FontSize',15)
    box off
    title('Variogram','Interpreter','Latex','FontSize',15)
   xlabel('Lag distance','Interpreter','Latex','FontSize', 15)
   ylabel('$\gamma (h)$','Interpreter','Latex','FontSize', 15)
   
   h=legend('Experimetal variogram','Fitted variogram');
set(h,'Interpreter','Latex','Location','southeast','FontSize',15);

end
%subplot(nplot,nplot,ndir+1);
%[nmodel,p]=size(model);
%if p==2
 %   tt1=strvcat('Mod. iso. : type, a, c');
%elseif deuxd==1
   % tt1=strvcat('Mod. aniso. 2D : type, ax ay, rot (trigo), c');
%else
 %   tt1=strvcat('Mod. aniso. 3D : type, ax, ay, az, rotx, roty, rotz, c');
%end

%tt4=strvcat('1-pépite','2-exponentiel','3-gaussien','4-sphérique');

%tt3=' ';
%for i=1:nmodel
 %   tt3=strvcat(tt3,tt4(model(i,1),:));
%end

%if deuxd==1
 %   tt2=num2str([model(:,[1:3,7]),c],4);
    
%else
 %   tt2=num2str([model,c],4);
%end
%tt2=[tt2(:,1),ones(nmodel,1)*'  ',tt2(:,2:end)];
%[n,p]=size(tt2);
%for i=p:-1:5;
 %   if strcmp(tt2(:,[i-2:i]),char(ones(n,3)*' '));
  %      tt2(:,i)=[];
  %  end
%end
%tt2=strvcat(tt1,tt2,tt3);

%t=text(0.20,1.15,tt1,'units','normalized','tag','model','horizontalalignment','left','fontname','courier','fontsize',7)
%t=text(0.10,0.50,tt2,'units','normalized','tag','model','horizontalalignment','left','fontname','courier','fontsize',7);
%axis off;
%hold off



% Ici, on définit des uicontrol pour pouvoir manipuler les objets du graphe
% On peut l'enlever si cela cause problème


if isempty(findobj(gcf,'label','Affichage')) 
    t=uimenu(gcf,'label','Affichage');
    t1=uimenu(t,'label','nombre de paires','checked','on','tag','apaire',...
        'callback','apaire');
    t2=uimenu(t,'label','modèle','checked','on','tag','amodele',...
        'callback','amodele2');
end
%set(gcf,'color',[1 1 1]);
figure(gcf)


function x=poletocart(pole)
%
% fonction pour convertir des pôles en coordonnées cartésiennes (système main droite)
% pole: nx2 dir (azimut), pendage du pole angle positif vers le bas (retourne une coordonnée z négative)
deg2rad=pi/180;
pole=pole*deg2rad;

x(:,1)=cos(pole(:,2)).*sin(pole(:,1));
x(:,2)=cos(pole(:,2)).*cos(pole(:,1));
x(:,3)=-sin(pole(:,2));

