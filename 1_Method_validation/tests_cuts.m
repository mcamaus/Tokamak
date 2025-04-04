R=0.5; Rp=0.3;

Currentshape="D";

sigma=2;  delta=0.33; beta=asin(delta); k=1.7;
Jrhoz=@(rho,z) Dcurrent(rho,z,R,Rp,k,beta,sigma);


rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);
alpha=linspace(-pi,pi);
dists=sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
maxdist=max(dists);
Rp=maxdist;


xvec=linspace(0.1,R+5*0.15,100); x=xvec;
zvec=linspace(-5*0.15,5*0.15,100); z=zvec;


figure;
plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'-r','LineWidth',2)

axis equal
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
shading interp
view(2)
axis equal
colormap(jet)
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');

axis([0.1,R+5*0.15,-5*0.15,5*0.15])
%% desde dentro
% Cálculo de alpha($Psi$)

hold off
xp=0.4; zp=0.2;
hold on
scatter(xp,zp,10)

Psi=pi/3;



H=@(alpha) sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
%H=@(alpha) Rp;

alphavec=linspace(0,2*pi,100);
plot(R+H(alphavec).*cos(atan2(zf(alphavec),rhof(alphavec)-R)),H(alphavec).*sin(atan2(zf(alphavec),rhof(alphavec)-R)))
axis equal


plot(rhof(alphavec),zf(alphavec),'--r')
axis equal

minfun=@(alpha) cos(Psi)*(H(alpha)*sin(alpha)-zp)-sin(Psi)*(H(alpha)*cos(alpha)+(R-xp));

alphasol=fminsearch(minfun,Rp);

Deltasupsol=((R-xp)*cos(Psi)-zp*sin(Psi))+sqrt(((R-xp)*cos(Psi)-zp*sin(Psi))^2+(H(alphasol)^2-zp^2-(R-xp)^2));

plot([xp,xp+Deltasupsol*cos(Psi)],[zp,zp+Deltasupsol*sin(Psi)],'-b')
%% 

sigma=2;  delta=0.33; beta=asin(delta); k=1.7;

rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);

theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));

% figure
% plot(wrapTo2Pi(thetavec),alphavec)


% hold on
% plot(thetavec,alphareg(thetavec),'x')
%Htheta=@(theta) Rp;
Htheta=@(theta) H(alphareg(theta));

plot(R+Htheta(alphavec).*cos(alphavec),Htheta(alphavec).*sin(alphavec))
axis equal


xp=0.4; zp=0.2;
hold on
scatter(xp,zp,10)
Psi=1*pi/6;


minfun=@(alpha) (cos(Psi)*(Htheta(alpha)*sin(alpha)-zp)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-xp)));
alphap=Psi+asin((sin(Psi)*(R-xp)+zp*cos(Psi))/Rp);
alphasol=fzero(minfun,alphap);
Deltasupsol=((R-xp)*cos(Psi)-zp*sin(Psi))+sqrt(((R-xp)*cos(Psi)-zp*sin(Psi))^2+(Htheta(alphasol)^2-zp^2-(R-xp)^2));
plot([xp,xp+Deltasupsol*cos(Psi)],[zp,zp+Deltasupsol*sin(Psi)],'-b')


%% desde fuera, limites D
sigma=2;  delta=0.33; beta=asin(delta); k=1.7;

rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);

theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));

% figure
% plot(wrapTo2Pi(thetavec),alphavec)


% hold on
% plot(thetavec,alphareg(thetavec),'x')
%Htheta=@(theta) Rp;
Htheta=@(theta) H(alphareg(theta));
Params=struct('Htheta',Htheta,'R',R);
Htheta=Params.Htheta; R=Params.R;
i=1;
figure
while i==1

plot(R+Htheta(alphavec).*cos(alphavec),Htheta(alphavec).*sin(alphavec))
axis equal

%xp=1; zp=0.2;
[xp,zp]=ginput(1);
hold on
scatter(xp,zp,10)
%Psi=pi-pi/4;

[Psimin,Psimax]=LimsPsiExteriorOpt(xp,zp,Params);
Psivec=linspace(Psimin,Psimax,100); Psi=Psivec(randi([1,100]));



Psic=atan2(-zp,(R-xp));
Deltac=sqrt((R-xp)^2+zp^2);

plot([xp,xp+Deltac*cos(Psic)],[zp,zp+Deltac*sin(Psic)],'--k');

minfun=@(alpha) abs(cos(Psi)*(Htheta(alpha)*sin(alpha)-zp)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-xp)));
alphaprop1=atan2(zp,xp-R);
alphaprop2=pi+alphaprop1;

feval1=1;
while feval1>0.01
[alphasol1,feval1]=fminsearch(minfun,alphap1);
alphaprop1=alphaprop1+rand*pi/2; n1=n1+1;
end


alphasol2=alphasol1;
feval2=1;n2=0;
while abs(wrapToPi(alphasol2-alphasol1))<1e-1 || feval2>0.01
[alphasol2,feval2]=fminsearch(minfun,alphaprop2);
alphaprop2=alphaprop2+rand*pi; n2=n2+1;
end


plot([xp,R+Htheta(alphasol1)*cos(alphasol1)],[zp,Htheta(alphasol1)*sin(alphasol1)],'o','MarkerSize',10,'LineWidth',2)
plot([xp,R+Htheta(alphasol2)*cos(alphasol2)],[zp,Htheta(alphasol2)*sin(alphasol2)],'x','MarkerSize',10,'LineWidth',2)
Deltamin=sqrt((R+Htheta(alphasol1)*cos(alphasol1)-xp)^2+(Htheta(alphasol1)*sin(alphasol1)-zp)^2);
Deltamax=sqrt((R+Htheta(alphasol2)*cos(alphasol2)-xp)^2+(Htheta(alphasol2)*sin(alphasol2)-zp)^2);


%Deltamax=Deltac*cos(Psi-Psic)+sqrt(Htheta(alphasol2)^2-Deltac^2*sin(Psi-Psic));

plot([xp,xp+Deltamin*cos(Psi)],[zp,zp+Deltamin*sin(Psi)])
plot([xp,xp+Deltamax*cos(Psi)],[zp,zp+Deltamax*sin(Psi)])

pause(2)
hold off
end
% 
% 
%% búsqueda de psimin psimax Rafa


Psic=atan2(-zp,(R-xp));
Deltac=sqrt((R-xp)^2+zp^2);

parameters=struct('R',R,'xp',xp,'zp',zp,'Htheta',Htheta);

minfun=@(Psi) abs(minfun_psi(Psi,parameters));

[~,Hmin]=fminsearch(Htheta,0);
[~,Hmax]=fminsearch(@(theta)-Htheta(theta),0);
Hmax=-Hmax;

Psiminp1=Psic-asin(Rp/Deltac);
Psiminp2=Psic-asin(Hmin/Deltac);

Psivec=linspace(Psic, Psiminp1,40);
difalpha=zeros(size(Psivec));

% for i=length(Psivec):-1:1
% difalpha(i)=minfun(Psivec(i));
% fprintf('Psivec= %5.4f, difalpha= %5.4f \n',Psivec(i),difalpha(i))
% %pause();
% end

alphavec=linspace(0,2*pi,100);
plot(R+Htheta(alphavec).*cos(alphavec),Htheta(alphavec).*sin(alphavec))
hold on; axis equal;
plot([xp,xp+1*cos(Psic)],[zp,zp+1*sin(Psic)])
plot([xp,xp+1*cos(Psiminp1)],[zp,zp+1*sin(Psiminp1)])
plot([xp,xp+1*cos(Psiminp2)],[zp,zp+1*sin(Psiminp2)])
axis([0.1,R+5*0.15,-5*0.15,5*0.15])

for i=length(Psivec):-1:1
    tic;
difalpha(i)=minfun(Psivec(i));
toc;
fprintf('Psivec= %5.4f, difalpha= %5.4f \n',Psivec(i),difalpha(i))
%pause();
end
figure
plot(Psivec,difalpha,'-x')

difalphaReg=@(Psi) pchip(Psivec(difalpha<100),difalpha(difalpha<100),Psi);

hold on
plot(Psivec,difalphaReg(Psivec),'-')

Psimin=fzero(difalphaReg,Psic);

plot([xp,xp+1*cos(Psimin)],[zp,zp+1*sin(Psimin)])
shg



axis([min(Psivec),max(Psivec),0,1])



%Psimin=fminsearch(minfun,Psiminp1);


%% busqueda de psimin psimax en función de las rectas tangentes
figure
i=1;
%while i<10 
alphavec=linspace(0,2*pi,100);
plot(R+Htheta(alphavec).*cos(alphavec),Htheta(alphavec).*sin(alphavec))
hold on; axis equal;
axis([0.1,R+5*0.15,-5*0.15,5*0.15])

xp=0.68; zp=0.45;
%[xp,zp]=ginput(1);

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

fprintf('Número intentos n1=%3i, n2=%3i; log10(feval1): %3.2f, log10(feval2): %3.2f \n',n1,n2,log10(feval1),log10(feval2))





plot(R+Htheta(alphasol1).*cos(alphasol1),Htheta(alphasol1).*sin(alphasol1),'x')
plot(R+Htheta(alphasol2).*cos(alphasol2),Htheta(alphasol2).*sin(alphasol2),'o')
%plot([xp,xp+1*cos(Psic)],[zp,zp+1*sin(Psic)])
plot([xp,R+Htheta(alphasol1).*cos(alphasol1)],[zp,Htheta(alphasol1).*sin(alphasol1)])
plot([xp,R+Htheta(alphasol2).*cos(alphasol2)],[zp,Htheta(alphasol2).*sin(alphasol2)])

plot([R+Htheta(alphasol2-deltaa).*cos(alphasol2-deltaa),R+Htheta(alphasol2+deltaa).*cos(alphasol2+deltaa)],...
    [Htheta(alphasol2-deltaa).*sin(alphasol2-deltaa),Htheta(alphasol2+deltaa).*sin(alphasol2+deltaa)])
Params=struct('Htheta',Htheta,'R',R);

deltax1=R+Htheta(alphasol1).*cos(alphasol1)-xp; deltaz1=Htheta(alphasol1).*sin(alphasol1)-zp;
deltax2=R+Htheta(alphasol2).*cos(alphasol2)-xp; deltaz2=Htheta(alphasol2).*sin(alphasol2)-zp;
Psisol1=wrapTo2Pi(atan2(deltaz1,deltax1)); Psisol2=wrapTo2Pi(atan2(deltaz2,deltax2));

Psimin1=min(Psisol1,Psisol2);
Psimax2=max(Psisol1,Psisol2);

[Psimin,Psimax]=LimsPsiExteriorOpt(xp,zp,Params);

title(sprintf('Psimin= %3.1f deg, Psimax=%3.1fdeg',Psimin*180/pi,Psimax*180/pi))
%pause()
hold off


%end


%%


function [Deltamin,Deltamax]=LimsDExtExact(xp,zp,Psi,Params)
Htheta=Params.Htheta; R=Params.R;
Psic=atan2(-zp,(R-xp));
Deltac=sqrt((R-xp)^2+zp^2);

minfun=@(alpha) abs(cos(Psi)*(Htheta(alpha)*sin(alpha)-zp)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-xp)));
alphaprop1=atan2(zp,xp-R);
alphaprop2=pi+alphaprop1;

feval1=1;
while feval1>0.01
[alphasol1,feval1]=fminsearch(minfun,alphap1);
alphaprop1=alphaprop1+rand*pi/2; n1=n1+1;
end


alphasol2=alphasol1;
feval2=1;n2=0;
while abs(wrapToPi(alphasol2-alphasol1))<1e-1 || feval2>0.01
[alphasol2,feval2]=fminsearch(minfun,alphaprop2);
alphaprop2=alphaprop2+rand*pi; n2=n2+1;
end

Delta1=sqrt((R+Htheta(alphasol1)*cos(alphasol1)-xp)^2+(Htheta(alphasol1)*sin(alphasol1)-zp)^2);
Delta2=sqrt((R+Htheta(alphasol2)*cos(alphasol2)-xp)^2+(Htheta(alphasol2)*sin(alphasol2)-zp)^2);

Deltamax=max(Delta1,Delta2);
Deltamin=min(Delta1,Delta2);


end

function [Psimin,Psimax]=LimsPsiExteriorOpt(xp,zp,Params)

Htheta=Params.Htheta; R=Params.R;

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
while abs(alphasol2-alphasol1)<1e-3 || feval2>0.01
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

%% funciones

function J=Dcurrent(rhov,zv,R,Rp,k,beta,sigma)

        rho=@(alpha,Rp) R+Rp/k*cos(alpha+beta*sin(alpha));
        z=@(alpha,Rp) Rp*sin(alpha);
        funz=@(x) norm([rhov-rho(x(1),x(2)),zv-z(x(1),x(2))]);

        alpha0=(atan2(zv,(rhov-R)));
        Rp0=sqrt((rhov-R)^2+zv^2);
        vec0=[alpha0;Rp0];
        vecmin=fminsearch(funz,vec0);

        Deltac=sqrt((rho(vecmin(1),Rp)-R)^2+z(vecmin(1),Rp)^2);
        Deltap=sqrt((rhov-R)^2+zv^2);

        Delta=Deltap/Deltac;

        J=exp(-sigma*(Deltap/Deltac).^2)*(Delta<1);



end