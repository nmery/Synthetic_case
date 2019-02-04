function [b,yc,e,sc,r2,s2b,rb]=regres(y,x,ord)
%
%function [b,yc,e,sc,r2,s2b,rb]=regres(y,x,ord)
%
% fonction effectuant la r‚gression de x sur y.
%
% Input: y: vecteur nx1, var. d‚pendante (… expliquer)
%        x: matrice nxp des var. ind‚pendantes
%        ord: 0, pas de constante dans le modŠle
%             1, constante dans le modŠle
%
% Output: b: vecteur (p+ord) des coefficients du modŠle, d'abord les p
%            coefficients des variables puis, si ord=1, la constante.
%         yc: vecteur nx1 avec les valeurs predites
%         e: vecteur nx1 des erreurs (yc-y)
%         sc: tableau Anova (1+ord)*3 x 3
%             ligne: SCT SCR SCE puis si ord=1 SCM SCTm SCRm
%             colonne: Somme des carr‚s, Carr‚ moyen, degr‚s de libert‚
%         r2: coefficient de correlation multiple 
%         s2b: matrice de variance-covariance des coefficients b
%         rb: matrice de corrélations entre les coefficients b
%         
[n,p]=size(x);
if ord==1,
   x=[x,ones(n,1)];
end;
    
b=x\y;				% coefficients de la regression
yc=x*b;				% valeurs predites
e=y-yc;				% erreurs
nl=(ord+1)*3;			% nombre de lignes dans sc selon le cas
sc=ones(nl,3)*nan;
sc(1,1)=y'*y;   sc(1,3)=n;	% SCT et d.l.
sc(2,1)=yc'*yc; sc(2,3)=p+ord;  % SCR et d.l.
sc(3,1)=e'*e;   sc(3,3)=sc(1,3)-sc(2,3);  % SCE et d.l.
if ord==1,
  sc(4,1)=n*mean(y)^2; sc(4,3)=1;         % SCM et d.l.
  sc(5,1)=sc(1,1)-sc(4,1) ; sc(5,3)=n-1;  % SCTm et d.l.	
  sc(6,1)=sc(2,1)-sc(4,1) ; sc(6,3)=p;    % SCRm et d.l.
end
sc(:,2)=sc(:,1)./max(sc(:,3),1);	% CM
if ord==1,
  r2=sc(6,1)/sc(5,1);		% R2	
else
  %  r2=sc(2,1)/sc(1,1);		      % ancienne façon de calculer R2 sans constante
  r2=1-sc(3,1)/(sc(1,1)-n*mean(y)^2);% nouvelle façon de calculer R2 sans constante
end
s2b=sc(3,2)*inv(x'*x);		% Variance-cov des coeff. b
t=sqrt(diag(s2b)); t=inv(diag(t)); 
rb=t*s2b*t;			 % correl. des coeff. b


% figure(999);
% plot(yc,y,'+b');hold on
% plot([min(y),max(y)],[min(y),max(y)],'-k','linewidth',2);
% title('Valeurs observées vs valeurs prédites','fontsize',14)
% xlabel('Prédites','fontsize',14)
% ylabel('Observées','fontsize',14)

hold off;%figure(999);


