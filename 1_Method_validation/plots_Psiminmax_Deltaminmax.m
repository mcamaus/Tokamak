
R=0.5; Rp=0.3;

 %rho=R+1.5*Rp; z=0*Rp;
 rho=R-1.1*Rp; z=-1.2*Rp;


%Currentshape="D";

sigma=2;  delta=0.33; beta=asin(delta); k=1.7; epsilon=Rp/k;

rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);

theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));
H=@(alpha) sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
Htheta=@(theta) H(alphareg(theta));
Rp1=max(Htheta(alphavec));
Htheta=@(theta) H(alphareg(theta))*Rp/Rp1;
%rho=R+0.1*Rp; z=0.1*Rp;
figure(1)

plot(R+Htheta(thetavec).*cos(thetavec),Htheta(thetavec).*sin(thetavec))
hold on
scatter(rho,z,10)
axis equal


Params=struct('Htheta',Htheta,'R',R,'Rp',Rp);

[Psimin,Psimax,Psic,Dc]=LimsPsiExteriorExact(rho,z,Params);


Psivec=linspace(Psimin,Psimax,100);
Dminvec1=zeros(size(Psivec));
Dmaxvec1=zeros(size(Psivec));
for i=1:length(Psivec)
    [Dmin,Dmax]=LimsDExtExact(rho,z,Psivec(i),Psimin,Psimax,Params);
    Dmaxvec1(i)=Dmax;
    Dminvec1(i)=Dmin;


plot([rho,rho+Dmin*cos(Psivec(i))],[z,z+Dmin*sin(Psivec(i))])
plot([rho,rho+Dmax*cos(Psivec(i))],[z,z+Dmax*sin(Psivec(i))])

end

figure(2)
cmap = colororder();
plot(Psivec,Dminvec1,Psivec,Dmaxvec1,'-','color', cmap(1,:)); hold on

%pause()


%
% Currenshape="invD"

delta=-0.3; beta=asin(delta); k=1.7; epsilon=Rp/k;
rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);

theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));
H=@(alpha) sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
Htheta=@(theta) H(alphareg(theta));
Rp1=max(Htheta(alphavec));
Htheta=@(theta) H(alphareg(theta))*Rp/Rp1;
%rho=R+0.1*Rp; z=0.1*Rp;
figure(1)

plot(R+Htheta(thetavec).*cos(thetavec),Htheta(thetavec).*sin(thetavec))
hold on
scatter(rho,z,10)
axis equal


Params=struct('Htheta',Htheta,'R',R,'Rp',Rp);

[Psimin2,Psimax2,Psic,Dc]=LimsPsiExteriorExact(rho,z,Params);


Psivec2=linspace(Psimin2,Psimax2,100);
Dminvec2=zeros(size(Psivec2));
Dmaxvec2=zeros(size(Psivec2));
for i=1:length(Psivec2)
    [Dmin,Dmax]=LimsDExtExact(rho,z,Psivec2(i),Psimin2,Psimax2,Params);
    Dmaxvec2(i)=Dmax;
    Dminvec2(i)=Dmin;


%plot([rho,rho+Dmin*cos(Psivec(i))],[z,z+Dmin*sin(Psivec(i))])
%plot([rho,rho+Dmax*cos(Psivec(i))],[z,z+Dmax*sin(Psivec(i))])

end

figure(2)




% Currenshape="Circle"

delta=-0.3; beta=asin(delta); k=1.7;
rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);

theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));
H=@(alpha) sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
Htheta=@(theta) Rp;
%Rp=max(Htheta(alphavec));
%rho=R+0.1*Rp; z=0.1*Rp;
figure(1)

plot(R+Htheta(thetavec).*cos(thetavec),Htheta(thetavec).*sin(thetavec))
hold on
scatter(rho,z,10)
axis equal


Params=struct('Htheta',Htheta,'R',R,'Rp',Rp);

[Psimin3,Psimax3,Psic,Dc]=LimsPsiExteriorExact(rho,z,Params);


Psivec3=linspace(Psimin3,Psimax3,100);
Dminvec3=zeros(size(Psivec3));
Dmaxvec3=zeros(size(Psivec3));
for i=1:length(Psivec3)
    [Dmin,Dmax]=LimsDExtExact(rho,z,Psivec3(i),Psimin3,Psimax3,Params);
    Dmaxvec3(i)=Dmax;
    Dminvec3(i)=Dmin;


%plot([rho,rho+Dmin*cos(Psivec(i))],[z,z+Dmin*sin(Psivec(i))])
%plot([rho,rho+Dmax*cos(Psivec(i))],[z,z+Dmax*sin(Psivec(i))])

end
%%
figure(2)
h1=plot([Psivec,flip(Psivec)],[Dminvec1,flip(Dmaxvec1)],'-','color', cmap(1,:),'LineWidth',3);
hold on
h2=plot([Psivec2,flip(Psivec2)],[Dminvec2,flip(Dmaxvec2)],'-','color', cmap(2,:),'LineWidth',3);
plot([Psivec3,flip(Psivec3)],[Dminvec3,flip(Dmaxvec3)],'-','color', cmap(3,:),'LineWidth',3);

plot([Psivec(end),Psivec(1)],[Dminvec1(end),Dmaxvec1(1)],'--','color', cmap(1,:),'HandleVisibility','off','LineWidth',3)
plot([Psivec2(end),Psivec2(1)],[Dminvec2(end),Dmaxvec2(1)],'--','color', cmap(2,:),'HandleVisibility','off','LineWidth',3)
plot([Psivec3(end),Psivec3(1)],[Dminvec3(end),Dmaxvec3(1)],'--','color', cmap(3,:),'HandleVisibility','off','LineWidth',3)
legend('PT D-shape','NT D-shape','Circle','Location','Northwest')
axes=gca;
axes.FontSize=20;
axes.FontName="Times New Roman";
axes.XLim=[0,1.6];
axes.YLim=[0.1,1.1];
xlabel('\Psi (rad)')
ylabel('\Delta_{max}, \Delta_{min} (m)')
%legend('PT D-shape','NT D-shape','Circle','Location','Northwest')

%%

%%
figure(3)

plot(Psivec,Dmaxvec1,Psivec,Dmaxvec2,Psivec,Dmaxvec3,'Linewidth',3)
legend('D shape','Inverted D shape','Circle','Location','best')

axes=gca;
axes.FontSize=20;
axes.FontName="Times New Roman";
axes.XLim=[0,2*pi];
xlabel('\Psi (rad)')
ylabel('\Delta_{sup} (m)')

%%

function [Deltamin,Deltamax]=LimsDExtExact(xp,zp,Psi,Psimin,Psimax,Params)
Htheta=Params.Htheta; R=Params.R;
Psic=atan2(-zp,(R-xp));
Deltac=sqrt((R-xp)^2+zp^2);

minfun=@(alpha) abs(cos(Psi)*(Htheta(alpha)*sin(alpha)-zp)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-xp)));
alphaprop1=atan2(zp,xp-R);
alphaprop2=pi+alphaprop1;

feval1=1;n1=0;
while feval1>0.01
[alphasol1,feval1]=fminsearch(minfun,alphaprop1);
alphaprop1=alphaprop1+rand*pi/2; n1=n1+1;
if n1>20
    warning('n1 LimdsDExtExact >20')
end
end


alphasol2=alphasol1;
feval2=1;n2=0;
while abs(wrapToPi(alphasol2-alphasol1))<1e-2 || feval2>0.01
[alphasol2,feval2]=fminsearch(minfun,alphaprop2);
alphaprop2=alphaprop2+rand*pi; n2=n2+1;
if abs(wrapToPi(Psimin-Psi))<0.01 || abs(wrapToPi(Psimax-Psi))<0.01
    break
elseif n2>20
    warning('n2 LimdsDExtExact >20')
end
end

Delta1=sqrt((R+Htheta(alphasol1)*cos(alphasol1)-xp)^2+(Htheta(alphasol1)*sin(alphasol1)-zp)^2);
Delta2=sqrt((R+Htheta(alphasol2)*cos(alphasol2)-xp)^2+(Htheta(alphasol2)*sin(alphasol2)-zp)^2);

Deltamax=max(Delta1,Delta2);
Deltamin=min(Delta1,Delta2);


end

function [Psimin,Psimax,Psic,Dc]=LimsPsiExteriorExact(xp,zp,Params)

Htheta=Params.Htheta; R=Params.R;

Psic=atan2(-zp,(R-xp));
Dc=sqrt((R-xp)^2+zp^2);
deltaa=0.001;
%xp; zp;
drx=@(alpha) (Htheta(alpha+deltaa)*cos(alpha+deltaa)-Htheta(alpha-deltaa)*cos(alpha-deltaa))/(2*deltaa);
drz=@(alpha) (Htheta(alpha+deltaa)*sin(alpha+deltaa)-Htheta(alpha-deltaa)*sin(alpha-deltaa))/(2*deltaa);

ecu=@(alpha)  abs((xp-(R+Htheta(alpha)*cos(alpha)))/drx(alpha)-(zp-Htheta(alpha)*sin(alpha))/drz(alpha));



alphaprop1=atan2(zp,xp-R)+pi/3;
alphaprop2=atan2(zp,xp-R)-pi/3;

feval1=1;n1=0;
while feval1>0.01
[alphasol1,feval1]=fminsearch(ecu,alphaprop1);
alphaprop1=alphaprop1+rand*pi/4; n1=n1+1;
end

alphasol2=alphasol1;
feval2=1;n2=0;
while abs(wrapToPi(alphasol2-alphasol1))<1e-3 || feval2>0.01
[alphasol2,feval2]=fminsearch(ecu,alphaprop2);
alphaprop2=alphaprop2+rand*pi/4; n2=n2+1;
end

deltax1=R+Htheta(alphasol1).*cos(alphasol1)-xp; deltaz1=Htheta(alphasol1).*sin(alphasol1)-zp;
deltax2=R+Htheta(alphasol2).*cos(alphasol2)-xp; deltaz2=Htheta(alphasol2).*sin(alphasol2)-zp;
Psisol1=wrapTo2Pi(atan2(deltaz1,deltax1)); Psisol2=wrapTo2Pi(atan2(deltaz2,deltax2));

Psimin=min(Psisol1,Psisol2);
Psimax=max(Psisol1,Psisol2);

if abs(Psimax-Psimin)>pi
    Psimax=min(Psisol1,Psisol2);
    Psimin=max(Psisol1,Psisol2)-2*pi;

end

end

function [difalpha]=minfun_psi(Psi,parameters)

R=parameters.R;
xp=parameters.xp;
zp=parameters.zp;
Htheta=parameters.Htheta;


Psic=atan2(-zp,(R-xp));
Deltac=sqrt((R-xp)^2+zp^2);



minfun=@(alpha) (cos(Psi)*(Htheta(alpha)*sin(alpha)-zp)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-xp)));
alphap1=Psic+pi;
alphap2=Psic;

options=optimset('MaxIter',10);

[alphasol1,feval1,exitflag,output]=fzero(minfun,alphap1,options);
if exitflag~=1
    difalpha=100;
    return
end
%minfun2=@(alpha) minfun(alpha)+0.5*(abs(alpha-alphasol1)<0.5);
[alphasol2,feval2,exitflag,output]=fzero(minfun,alphap2,options);

if 0
    alphavec=linspace(0,2*pi,100);
    plot(R+Htheta(alphavec).*cos(alphavec),Htheta(alphavec).*sin(alphavec))
    hold on; axis equal;
    plot([xp,R+Htheta(alphasol1)*cos(alphasol1)],[zp,Htheta(alphasol1)*sin(alphasol1)])
    plot([xp,R+Htheta(alphasol2)*cos(alphasol2)],[zp,Htheta(alphasol2)*sin(alphasol2)])
    plot([R,R+Htheta(alphasol1)*cos(alphasol1)],[0,Htheta(alphasol1)*sin(alphasol1)])
    plot([R,R+Htheta(alphasol2)*cos(alphasol2)],[0,Htheta(alphasol2)*sin(alphasol2)])
    plot([xp,xp+10*cos(Psi)],[zp,zp+10*sin(Psi)])
   
    axis([0.1,R+5*0.15,-5*0.15,5*0.15])
    shg
    hold off
end

if max(feval1,feval2)>0.01
alphasol1=0; alphasol2=1;
end

difalpha=wrapToPi(alphasol1-alphasol2);
end

function [Dmin,Dmax]=LimsDIntExact(rho,z,Psi,Params)
R=Params.R; Rp=Params.Rp; Htheta=Params.Htheta;


Dmin=0;
minfun=@(alpha) (cos(Psi)*(Htheta(alpha)*sin(alpha)-z)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-rho)));
alphap=Psi+asin((sin(Psi)*(R-rho)+z*cos(Psi))/Rp);
alphasol=fzero(minfun,alphap);
Dmax=((R-rho)*cos(Psi)-z*sin(Psi))+sqrt(((R-rho)*cos(Psi)-z*sin(Psi))^2+(Htheta(alphasol)^2-z^2-(R-rho)^2));

end
%% funciones