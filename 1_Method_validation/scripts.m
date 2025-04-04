
%% Fig 4(a)

%R=0.5; Rp=0.15;
R=0.5; Rp=0.15;
Currentshape="Constant";
%Currentshape="Gaussian"; 
%Currentshape="Kidney"; 
%Currentshape="2SignKid"; 


if Currentshape=="Constant"
    a=1; b=1;
Jrhoz=@(rho,z) 1*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
Htheta=@(rho,z) Rp;
elseif Currentshape=="Gaussian"
sigma=0.75; a=1; b=1;
Jrhoz=@(rho,z) exp(-(rho-R).^2/(a*sigma*Rp)^2-z.^2/(b*sigma*Rp)^2).*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Kidney"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid(rho,z,a,b,sigma,R,Rp);
elseif Currentshape=="2SignKid"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid2sign(rho,z,a,b,sigma,R,Rp);
end

vecNp=(5:5:40);
Params=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',5,'DeltaDif',1e-6, 'test',0,'Htheta',Htheta,'ExactLims',0);


BNpx=zeros(size(vecNp)); BNpz=BNpx; ANp=BNpx; BANpz=BNpz; BANpx=BNpx; ADirNp=ANp; ADE=ANp; ADExy=ADE;

xp=R+0.1*Rp; zp=0.1*Rp;


for i=1:length(vecNp)

Params.MRW=vecNp(i);
%[BNpx(i),BNpz(i)]=FieldB(xp,zp,Params);
%[BANpx(i),BANpz(i)]=BfromA(xp,zp,Params);
ANp(i)=VectorPot(xp,zp,Params);
ADirNp(i)=VectorPotDirMRW(xp,zp,Params);
ADE(i)=VectorPotDE6(xp,zp,Params,i-1);
ADExy(i)=VectorPotDE(xp,zp,Params,i-1);
end

%M=6;
%ADE(end)=VectorPotDE2(xp,zp,Params,M);

names={'Original';'MRW Direct'; 'Doub. Exp. in polar';'Doub. Exp. in xy' };
sols=[ANp(end); ADirNp(end); ADE(end); ADExy(end)];
table(names,sols)
 %return


% subplot(1,2,1)
% title('Direct B')
% semilogy(vecNp,abs(BNpx-BNpx(end))/abs(BNpx(end)),vecNp,abs(BNpz-BNpz(end))/abs(BNpz(end)),vecNp,abs(ANp-ANp(end))/abs(ANp(end)))
% subplot(1,2,2)
% title('B from A')
% semilogy(vecNp,abs(BANpx-BANpx(end))/abs(BANpx(end)),vecNp,abs(BANpz-BANpz(end))/abs(BANpz(end)),vecNp,abs(ANp-ANp(end))/abs(ANp(end)))
% legend('B_x (R,R_p/2)','B_z (R,R_p/2)','A_\phi (R,R_p/2)')

figure
t = tiledlayout(1,1);
ax1 = axes(t);
h=plot(ax1,vecNp,log10(abs(ANp-ANp(end))/abs(ANp(end))),'b-','LineWidth',3); hold on
%h2=plot(ax1,vecNp,log10(abs(ADirNp-ADirNp(end))/abs(ADirNp(end))),'r-','LineWidth',3);
ax2=axes(t);
h3=plot(ax2,0:(length(vecNp)-1),log10(abs(ADE-ADE(end))/abs(ADE(end))),'r-','LineWidth',3);
ax2.XAxisLocation = 'top';
%ax1.Color = 'none';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlabel(ax1,'N (MRW)'); ylabel(ax1,'$\log_{10}|(A_\varphi -A_\varphi^\mathrm{ref})/A_\varphi^\mathrm{ref}|$','Interpreter','latex')
xlabel(ax2,'M (DE)')
legend([h,h3],'$A_\varphi$ (MRW)','$A_\varphi$ (DE)','Interpreter','latex','Location','Southwest')
ax1.FontSize=25;
ax2.FontSize=25;
%ax.YTick=-15:2:0;
ax1.FontName='Times New Roman';
ax2.FontName='Times New Roman';
%title('R=0.5 m, R_p=0.15 m')
axis(ax1,[5,35,-16,0])
axis(ax2,[0,6,-16,0])
box on
shg

%% Fig. 4(b)

Params.R=0.5; Params.Rp=0.15;
vecNp=(5:5:40);
BNpx=zeros(size(vecNp)); BNpz=BNpx; ANp=BNpx; BANpz=BNpz; BANpx=BNpx; AdirNp=ANp; ADE=ANp;

xp=R+0.1*Rp; zp=0.1*Rp;
for i=1:length(vecNp)

Params.MRW=vecNp(i);
[BNpx(i),BNpz(i)]=FieldB(xp,zp,Params);
[BANpx(i),BANpz(i)]=BfromA(xp,zp,Params);
[BNpxDir(i),BNpzDir(i)]=FieldBDirMRW(xp,zp,Params);
[BNpxDE(i),BNpzDE(i)]=FieldBDE(xp,zp,Params,i-1);
end

%M=6;
%ADE(end)=VectorPotDE(xp,zp,Params,M);

 

names={'Original';'MRW Direct'; 'Doub. Exp. in distance'};
solx=[BNpx(end); BNpxDir(end); BNpxDE(end)];
solz=[BNpz(end); BNpzDir(end); BNpzDE(end)];
table(names,solx,solz)

figure
t = tiledlayout(1,1);
ax1 = axes(t);
h=plot(ax1,vecNp,log10(abs(BNpx-BNpx(end))/abs(BNpx(end))),'b-','LineWidth',3); hold on
%h2=plot(ax1,vecNp,log10(abs(BNpxDir-BNpxDir(end))/abs(BNpxDir(end))),'r-','LineWidth',3);
h3=plot(ax1,vecNp,log10(abs(BNpz-BNpz(end))/abs(BNpz(end))),'b--','LineWidth',3); hold on
%h4=plot(ax1,vecNp,log10(abs(BNpzDir-BNpzDir(end))/abs(BNpzDir(end))),'r--','LineWidth',3);
ax2=axes(t);
h5=plot(ax2,0:7,log10(abs(BNpxDE-BNpxDE(end))/abs(BNpxDE(end))),'r-','LineWidth',3); hold on
h6=plot(ax2,0:7,log10(abs(BNpzDE-BNpzDE(end))/abs(BNpzDE(end))),'r--','LineWidth',3);
ax2.XAxisLocation = 'top';
%ax1.Color = 'none';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlabel(ax1,'N (MRW)'); ylabel(ax1,'$\log_{10}|(B_{\rho,z} -B_{\rho,z}^\mathrm{ref})/B_{\rho,z}^\mathrm{ref}|$','Interpreter','latex')
xlabel(ax2,'M (DE)')
legend([h,h3,h5,h6],'$B_\rho$ (MRW)','$B_z $ (MRW)','$B_\rho$ (DE)',' $B_z$ (DE)','Interpreter','latex','Location','Southwest')

ax1.FontSize=25;
ax2.FontSize=25;
%ax.YTick=-15:2:0;
ax1.FontName='Times New Roman';
ax2.FontName='Times New Roman';
axis(ax1,[5,35,-16,0])
axis(ax2,[0,6,-16,0])
%title('R=0.5 m, R_p=0.15 m')
%axis([5,35,-15,-3])
box on
shg
%% Fig 4(c)

%R=0.5; Rp=0.15;
R=0.5; Rp=0.15;
Currentshape="Constant";
%Currentshape="Gaussian"; 
%Currentshape="Kidney"; 
%Currentshape="2SignKid"; 


if Currentshape=="Constant"
    a=1; b=1;
Jrhoz=@(rho,z) 1*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
Htheta=@(rho,z) Rp
elseif Currentshape=="Gaussian"
sigma=0.75; a=1; b=1;
Jrhoz=@(rho,z) exp(-(rho-R).^2/(a*sigma*Rp)^2-z.^2/(b*sigma*Rp)^2).*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Kidney"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid(rho,z,a,b,sigma,R,Rp);
elseif Currentshape=="2SignKid"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid2sign(rho,z,a,b,sigma,R,Rp);
end




vecNp=(5:5:40);
Params=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',5,'DeltaDif',1e-6, 'test',0,'Htheta',Htheta);


BNpx=zeros(size(vecNp)); BNpz=BNpx; ANp=BNpx; BANpz=BNpz; BANpx=BNpx; ADirNp=ANp; ADE=ANp; ADExy=ADE;

xp=R+0.9*Rp; zp=0.9*Rp;


for i=1:length(vecNp)

Params.MRW=vecNp(i);
%[BNpx(i),BNpz(i)]=FieldB(xp,zp,Params);
%[BANpx(i),BANpz(i)]=BfromA(xp,zp,Params);
ANp(i)=VectorPot(xp,zp,Params);
ADirNp(i)=VectorPotDirMRW(xp,zp,Params);
ADE(i)=VectorPotDE6(xp,zp,Params,i-1);
ADExy(i)=VectorPotDE(xp,zp,Params,i-1);
end

%M=6;
%ADE(end)=VectorPotDE2(xp,zp,Params,M);

names={'Original';'MRW Direct'; 'Doub. Exp. in polar';'Doub. Exp. in xy' };
sols=[ANp(end); ADirNp(end); ADE(end); ADExy(end)];
table(names,sols)
 %return


% subplot(1,2,1)
% title('Direct B')
% semilogy(vecNp,abs(BNpx-BNpx(end))/abs(BNpx(end)),vecNp,abs(BNpz-BNpz(end))/abs(BNpz(end)),vecNp,abs(ANp-ANp(end))/abs(ANp(end)))
% subplot(1,2,2)
% title('B from A')
% semilogy(vecNp,abs(BANpx-BANpx(end))/abs(BANpx(end)),vecNp,abs(BANpz-BANpz(end))/abs(BANpz(end)),vecNp,abs(ANp-ANp(end))/abs(ANp(end)))
% legend('B_x (R,R_p/2)','B_z (R,R_p/2)','A_\phi (R,R_p/2)')

figure
t = tiledlayout(1,1);
ax1 = axes(t);
h=plot(ax1,vecNp,log10(abs(ANp-ANp(end))/abs(ANp(end))),'b-','LineWidth',3); hold on
%h2=plot(ax1,vecNp,log10(abs(ADirNp-ADirNp(end))/abs(ADirNp(end))),'r-','LineWidth',3);
ax2=axes(t);
h3=plot(ax2,0:(length(vecNp)-1),log10(abs(ADE-ADE(end))/abs(ADE(end))),'r-','LineWidth',3);
ax2.XAxisLocation = 'top';
%ax1.Color = 'none';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlabel(ax1,'N (MRW)'); ylabel(ax1,'$\log_{10}|(A_\varphi -A_\varphi^\mathrm{ref})/A_\varphi^\mathrm{ref}|$','Interpreter','latex')
xlabel(ax2,'M (DE)')
legend([h,h3],'$A_\varphi$ (MRW)','$A_\varphi$ (DE)','Interpreter','latex','Location','Southwest')
ax1.FontSize=25;
ax2.FontSize=25;
%ax.YTick=-15:2:0;
ax1.FontName='Times New Roman';
ax2.FontName='Times New Roman';
%title('R=0.5 m, R_p=0.15 m')
axis(ax1,[5,35,-16,0])
axis(ax2,[0,6,-16,0])
box on
shg

%% Fig. 4(d)

Params.R=0.5; Params.Rp=0.15;
vecNp=(5:5:40);
BNpx=zeros(size(vecNp)); BNpz=BNpx; ANp=BNpx; BANpz=BNpz; BANpx=BNpx; AdirNp=ANp; ADE=ANp;

xp=R+0.9*Rp; zp=0.9*Rp;
for i=1:length(vecNp)

Params.MRW=vecNp(i);
[BNpx(i),BNpz(i)]=FieldB(xp,zp,Params);
[BANpx(i),BANpz(i)]=BfromA(xp,zp,Params);
[BNpxDir(i),BNpzDir(i)]=FieldBDirMRW(xp,zp,Params);
[BNpxDE(i),BNpzDE(i)]=FieldBDE(xp,zp,Params,i-1);
end

%M=6;
%ADE(end)=VectorPotDE(xp,zp,Params,M);

 

names={'Original';'MRW Direct'; 'Doub. Exp. in distance'};
solx=[BNpx(end); BNpxDir(end); BNpxDE(end)];
solz=[BNpz(end); BNpzDir(end); BNpzDE(end)];
table(names,solx,solz)

figure
t = tiledlayout(1,1);
ax1 = axes(t);
h=plot(ax1,vecNp,log10(abs(BNpx-BNpx(end))/abs(BNpx(end))),'b-','LineWidth',3); hold on
%h2=plot(ax1,vecNp,log10(abs(BNpxDir-BNpxDir(end))/abs(BNpxDir(end))),'r-','LineWidth',3);
h3=plot(ax1,vecNp,log10(abs(BNpz-BNpz(end))/abs(BNpz(end))),'b--','LineWidth',3); hold on
%h4=plot(ax1,vecNp,log10(abs(BNpzDir-BNpzDir(end))/abs(BNpzDir(end))),'r--','LineWidth',3);
ax2=axes(t);
h5=plot(ax2,0:7,log10(abs(BNpxDE-BNpxDE(end))/abs(BNpxDE(end))),'r-','LineWidth',3); hold on
h6=plot(ax2,0:7,log10(abs(BNpzDE-BNpzDE(end))/abs(BNpzDE(end))),'r--','LineWidth',3);
ax2.XAxisLocation = 'top';
%ax1.Color = 'none';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlabel(ax1,'N (MRW)'); ylabel(ax1,'$\log_{10}|(B_{\rho,z} -B_{\rho,z}^\mathrm{ref})/B_{\rho,z}^\mathrm{ref}|$','Interpreter','latex')
xlabel(ax2,'M (DE)')
legend([h,h3,h5,h6],'$B_\rho$ (MRW)','$B_z $ (MRW)','$B_\rho$ (DE)',' $B_z$ (DE)','Interpreter','latex','Location','Southwest')

ax1.FontSize=25;
ax2.FontSize=25;
%ax.YTick=-15:2:0;
ax1.FontName='Times New Roman';
ax2.FontName='Times New Roman';
axis(ax1,[5,35,-16,0])
axis(ax2,[0,6,-16,0])
%title('R=0.5 m, R_p=0.15 m')
%axis([5,35,-15,-3])
box on
shg

%% Fig 5ab

x=linspace(0.1,1.5,100);
z=linspace(-1.5,1.5,100);
tic;
[f,A1,BDir]=PlotfwB(VecPotfun,FieldBfun,x,z,Params);
toc;
tic;
[f,A2,BA]=PlotfwB(VecPotfun,BfromAfun,x,z,Params);
toc;

ErrorBx=log(abs(BDir.Bx-BA.Bx)./abs(BDir.Bx));
ErrorBz=log(abs(BDir.Bz-BA.Bz)./abs(BDir.Bz));
%% 5a
figure
surf(x,z,0*ErrorBx,ErrorBx); hold on;
view(2)
shading flat
colormap(viridis)
ax=gca;
ax.FontSize=20;
xlabel("\rho (m)")
ylabel("z (m)")
ax.FontName='Times New Roman';
%title('R=0.5 m, R_p=0.15 m')
h=colorbar;
h.Ticks=-28:3:-4;
ylabel(h,"$\log_{10}|(B_\rho^{der}-B_\rho^{int})/B_\rho^{int}|$", 'Interpreter','latex');
axis equal
axis([min(x),max(x),min(z),max(z)])
%% 5b
figure
surf(x,z,0*ErrorBx,ErrorBz); hold on;
view(2)
shading flat
colormap(viridis)
ax=gca;
ax.FontSize=20;
xlabel("\rho (m)")
ylabel("z (m)")
ax.FontName='Times New Roman';
%title('R=0.5 m, R_p=0.15 m')
h=colorbar;
h.Ticks=-28:3:-4;
ylabel(h,"$\log_{10}|(B_z^{der}-B_z^{int})/B_z^{int}|$", 'Interpreter','latex');
axis equal
axis([min(x),max(x),min(z),max(z)])

%% Fig 6ab

R=0.5; Rp=0.001;
Jrhoz=@(rho,z) 1/(pi*Rp^2)*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
Params1=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',30,'DeltaDif',1e-6, 'test',0);

x=linspace(0.35,0.65,300); z=linspace(-0.15,0.15,300);

Ae=zeros(length(x),length(z));
Bxe=Ae; Bze=Ae;
A=Ae; Bx=Ae; Bz=Ae;

vec1=1:length(x); vec2=1:length(z);
tic;
parfor i=vec1
    for j=vec2

A(j,i)=VectorPot(x(i),z(j),Params1);
[Bx(j,i),Bz(j,i)]=FieldB(x(i),z(j),Params1);

[Ae(j,i),Bxe(j,i),Bze(j,i)]=FieldLoop(x(i),z(j),Params1);


    end
end
toc;
%% 6a

figure
surf(x,z,0*A,A*1e6); hold on;
%h1=quiver(x,z,Bx,Bz,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(x,z,Bx,Bz);
set(h1,'Color','k')
axis equal
axis([0.35,0.65,-0.15,0.1501])
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)','Interpreter','latex');
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Tokamak Rp=1 mm")
legend(h1,'Magnetic field','Location','Northeast')

%% 6b

figure
surf(x,z,0*A,Ae*1e6); hold on;
%h1=quiver(x,z,Bx,Bz,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(x,z,Bxe,Bze);
set(h1,'Color','k')
axis equal
axis([0.35,0.65,-0.15,0.1501])
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)','Interpreter','latex');
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Tokamak Rp=1 mm")
legend(h1,'Magnetic field','Location','Northeast')

%% Fig 7ab

R=50; Rp=0.25;
Jrhoz=@(rho,z) 1/(pi*Rp^2)*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
Params1=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',30,'DeltaDif',1e-6, 'test',0);

x=linspace(49.3,50.7,300); z=linspace(-0.7,0.7,300);

Ae=zeros(length(x),length(z));
Bxe=Ae; Bze=Ae;
A=Ae; Bx=Ae; Bz=Ae;

vec1=1:length(x); vec2=1:length(z);
tic;
parfor i=vec1
    for j=vec2

A(j,i)=VectorPot(x(i),z(j),Params1);
[Bx(j,i),Bz(j,i)]=FieldB(x(i),z(j),Params1);

[Ae(j,i),Bxe(j,i),Bze(j,i)]=FieldCilinder(x(i),z(j),Params1);


    end
end

Ae=Ae/(pi*Rp^2); Bxe=Bxe/(pi*Rp^2); Bze=Bze/(pi*Rp^2);
toc;
%% Fig7a

figure
surf(x,z,0*A,A*1e6); hold on;
%h1=quiver(x,z,Bx,Bz,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(x,z,Bx,Bz);
set(h1,'Color','k')
axis equal
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)', 'Interpreter','latex');
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Tokamak R=50 m, Rp=0.25 m")
legend(h1,'Magnetic field','Location','Northeast')


%% Fig 7b

dif=max(max(A))-max(max(Ae));
figure
surf(x,z,0*Ae,(Ae+dif)*1e6); hold on;
%h1=quiver(x,z,Bxe,Bze,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(x,z,Bxe,Bze);
set(h1,'Color','k')
axis equal
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)', 'Interpreter','latex');
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Infinite cylinder")
legend(h1,'Magnetic field','Location','Northeast')


%% Fig 7c d

R=50; Rp=0.25;
Jrhoz=@(rho,z) 1/(pi*Rp^2)*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
Params1=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',30,'DeltaDif',1e-6, 'test',0);

x=linspace(49.3,50.7,300); z=0;

Ae=zeros(length(x),length(z));
Bxe=Ae; Bze=Ae;
A=Ae; Bx=Ae; Bz=Ae;

vec1=1:length(x); vec2=1:length(z);
tic;
parfor i=vec1
 
A(i)=VectorPot(x(i),z,Params1);
[Bx(i),Bz(i)]=FieldB(x(i),z,Params1);

[Ae(i),Bxe(i),Bze(i)]=FieldCilinder(x(i),z,Params1);
    
end
toc;

%% Fig 7c
figure
dif=max(A)-max(Ae);
plot(x,A*1e6,'LineWidth',2); hold on
plot(x(1:8:end),1e6*(Ae(1:8:end)+dif),'+','LineWidth',3,'MarkerSize',15);
xlabel('\rho (m)');
ylabel('$A_\varphi$ (Tm)','Interpreter','latex');
legend('Torus', 'Inf. cylinder')
axis([0.9999*min(x),1.0001*max(x),min(A*1e6),1.01*max(A*1e6)])
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=25;
%% Fig 7d
figure

plot(x,Bz*1e6,'LineWidth',2); hold on
plot(x(1:8:end),Bze(1:8:end)*1e6,'+','MarkerSize',15,'LineWidth',2);
xlabel('\rho (m)');
ylabel('B_z (T)');
legend('Torus', 'Inf. cylinder')
axis([0.9999*min(x),1.0001*max(x),min(Bz)*1e6,1.01*max(Bz)*1e6])
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=25;

%% Fig 8a
R=0.5; Rp=0.3;

%R=0.5; Rp=0.001;

%Currentshape="Constant";
%Currentshape="Gaussian"; 
%Currentshape="Kidney"; 
%Currentshape="2SignKid"; 
Currentshape="D";

if Currentshape=="Constant"
    a=1; b=1;
Jrhoz=@(rho,z) 1*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Gaussian"
sigma=0.75;
a=1;
b=2;
Jrhoz=@(rho,z) exp(-(rho-R).^2/(a*sigma*Rp)^2-z.^2/(b*sigma*Rp)^2).*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Kidney"
sigma=0.1; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid(rho,z,a,b,sigma,R,Rp);%.*(rho>0);
elseif Currentshape=="2SignKid"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid2sign(rho,z,a,b,sigma,R,Rp);
elseif Currentshape=="D"
   sigma=2;  delta=0.33; beta=asin(delta); k=1.7;

 Jrhoz=@(rho,z) Dcurrent(rho,z,R,Rp,k,beta,sigma);



rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);
alpha=linspace(-pi,pi);


theta=@(alpha) atan2(zf(alpha),rhof(alpha)-R);

alphavec=linspace(0,2*pi,100);
thetavec=theta(alphavec);
alphareg=@(theta) pchip(wrapTo2Pi(thetavec),alphavec,wrapTo2Pi(theta));
H=@(alpha) sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
Htheta=@(theta) H(alphareg(theta));

dists=sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
maxdist=max(dists);
Rp=maxdist;

Params=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',15,'DeltaDif',1e-6, 'test',0,'Htheta',Htheta);
end


Params.ExactLims=1;
Params.MRWn=load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW)).nodes;
Params.MRWw=load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW)).weights;

Intfun=@(r) IntegrandH(r,0,R,0,Params);


VecPotfun=@(rho,z) VectorPot(rho,z,Params);
BfromAfun=@(rho,z) BfromA(rho,z,Params);
FieldBfun=@(rho,z) FieldB(rho,z,Params);

xvec=linspace(0.1,R+5*0.15,50); x=xvec;
zvec=linspace(-5*0.15,5*0.15,50); z=zvec;


[xmat,zmat]=ndgrid(xvec,zvec);
Jrhom=zeros(size(xmat));
for i=1:length(xvec)
    for j=1:length(zvec)

        Jrhom(i,j)=Jrhoz(xvec(i),zvec(j));



    end
end
hold on

surf(xmat,zmat,Jrhom)
hold on
%theta=linspace(0,2*pi);
%rho=Rp*(1-0.5*sin(theta/2).^15);
%plot3(R+rho.*cos(theta)*a,rho.*sin(theta)*b,1*ones(size(theta)),'--r','Linewidth',2);
%plot3(R+Rp.*cos(theta),Rp.*sin(theta),1*ones(size(theta)),'--r','Linewidth',2);

plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'--r','LineWidth',3)




axis equal
colormap(jet)
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
h=colorbar;
ylabel(h,'$J_\varphi$ (MA/$\mathrm{m}^2$)','Interpreter','latex');
caxis([0,1])
axis([0.1,R+5*0.15,-5*0.15,5*0.15])

%% Fig 8b

vec1=1:length(xvec); vec2=1:length(zvec);
tic;
A=zeros(length(xvec),length(zvec));
Bx=A; Bz=A;
problematic=[0,0];

fid=fopen('log.txt','w');
for i=vec1
    for j=vec2

%fprintf('x= %3.2f, z=%3.2f \n', xvec(i),zvec(j))
tic;
A(j,i)=VectorPot(xvec(i),zvec(j),Params);
[Bx(j,i),Bz(j,i)]=FieldB(xvec(i),zvec(j),Params);
toc;

%fprintf('log10s A: %3.2f, Bx= %3.2f, Bz= %3.2f \n', log10(A(j,i)),log10(Bx(j,i)),log10(Bz(j,i)));
    end
end

%%
%It=integral2(Jrhoz,R-Rp,R+Rp,-Rp,Rp);
figure
surf(xvec,zvec,0*A,A*1e6); hold on;
%h1=quiver(x,z,Bx,Bz,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(xvec,zvec,Bx,Bz);
set(h1,'Color','k')
axis equal
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)', 'Interpreter','latex');
%plot3(R+Rp*a*cos(t),Rp*b*sin(t),1*ones(size(t)),'r--','LineWidth',2);
plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'--r','LineWidth',3)
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Tokamak R=50 m, Rp=0.25 m")




%% Fig 9a

R=0.5; Rp=0.42;


%R=0.5; Rp=0.001;

%Currentshape="Constant";
%Currentshape="Gaussian"; 
%Currentshape="Kidney"; 
%Currentshape="2SignKid"; 
Currentshape="invD";

if Currentshape=="Constant"
    a=1; b=1;
Jrhoz=@(rho,z) 1*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Gaussian"
sigma=0.75;
a=1;
b=2;
Jrhoz=@(rho,z) exp(-(rho-R).^2/(a*sigma*Rp)^2-z.^2/(b*sigma*Rp)^2).*(((rho-R).^2/a^2+(z).^2/b^2)<=Rp^2);
elseif Currentshape=="Kidney"
sigma=0.1; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid(rho,z,a,b,sigma,R,Rp);%.*(rho>0);
elseif Currentshape=="2SignKid"
sigma=0.2; a=0.75; b=2; 
Jrhoz=@(rho,z) Jkid2sign(rho,z,a,b,sigma,R,Rp);
elseif Currentshape=="D"
   sigma=2;  delta=0.33; beta=asin(delta); k=1.7;

 Jrhoz=@(rho,z) Dcurrent(rho,z,R,Rp,k,beta,sigma);

elseif Currentshape=="invD"
load("nt_data.mat")

j=j*1e-6;

Jrhoz=@(rhop,zp) interp2(r,z,j,rhop,zp);



end




Params=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',15,'DeltaDif',1e-6, 'test',0);


Params.MRWn=load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW)).nodes;
Params.MRWw=load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW)).weights;

Intfun=@(r) IntegrandH(r,0,R,0,Params);


VecPotfun=@(rho,z) VectorPot(rho,z,Params);
BfromAfun=@(rho,z) BfromA(rho,z,Params);
FieldBfun=@(rho,z) FieldB(rho,z,Params);

xvec=linspace(0.1,R+5*0.15,100); x=xvec;
zvec=linspace(-5*0.15,5*0.15,100); z=zvec;


[xmat,zmat]=ndgrid(xvec,zvec);
Jrhom=zeros(size(xmat));
for i=1:length(xvec)
    for j=1:length(zvec)

        Jrhom(i,j)=Jrhoz(xvec(i),zvec(j));



    end
end
hold on
%% Plot Fig 9a
surf(xmat,zmat,Jrhom)
hold on
%theta=linspace(0,2*pi);
%rho=Rp*(1-0.5*sin(theta/2).^15);
%plot3(R+rho.*cos(theta)*a,rho.*sin(theta)*b,1*ones(size(theta)),'--r','Linewidth',2);
%plot3(R+Rp.*cos(theta),Rp.*sin(theta),1*ones(size(theta)),'--r','Linewidth',2);

%plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'--r','LineWidth',3)




axis equal
colormap(jet)
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
h=colorbar;
ylabel(h,'$J_\varphi$ (MA/$\mathrm{m}^2$)','Interpreter','latex');
%caxis([0,1])
axis([0.1,1.01,-5*0.15,5*0.15])
shg

%% Fig 9b,d
vec1=1:length(xvec); vec2=1:length(zvec);
tic;
A=zeros(length(xvec),length(zvec));
Bx=A; Bz=A;
for i=vec1
    for j=vec2
tic;
A(j,i)=VectorPot(xvec(i),zvec(j),Params);
toc;
[Bx(j,i),Bz(j,i)]=FieldB(xvec(i),zvec(j),Params);

    end
end
toc

%% Plot 9b
%It=integral2(Jrhoz,R-Rp,R+Rp,-Rp,Rp);
figure
surf(xvec,zvec,0*A,A*1e6); hold on;
%h1=quiver(x,z,Bx,Bz,'r');
%streamline(x,z,B.Bx,B.Bz,linspace(min(x),R+Rp,10),0*linspace(1,2,10))
h1=streamslice(xvec,zvec,Bx,Bz);
set(h1,'Color','k')
axis equal
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
colormap(jet)
h=colorbar;
ylabel(h,'$A_\varphi$ (Tm)', 'Interpreter','latex');
%plot3(R+Rp*a*cos(t),Rp*b*sin(t),1*ones(size(t)),'r--','LineWidth',2);
%plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'--r','LineWidth',3)
shading interp
view(2)
xlabel('\rho (m)'); ylabel('z(m)');
%title("Tokamak R=50 m, Rp=0.25 m")


%% Fig 9c

load('nt_data.mat')
figure
surf(r,z,Btheta_plasma)
shg
%figure; surf(xvec,zvec,sqrt(Bx.^2+Bz.^2))
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
h=colorbar;
ylabel(h,'$|B|$ (T)','Interpreter','latex');
shg
shading interp
view(2)
axis([0.1,1,-0.5,0.5])
colormap(jet)
xlabel('\rho (m)'); ylabel('z (m)')

%% Fig 9d

figure; surf(xvec,zvec,sqrt(Bx.^2+Bz.^2))
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
h=colorbar;
ylabel(h,'$|B|$ (T)','Interpreter','latex');
shg
shading interp
view(2)
axis([0.1,1,-0.5,0.5])
colormap(jet)
xlabel('\rho (m)'); ylabel('z (m)')


%% functions
function [A,Integrand,nodesr,nodesz]=VectorPotDirTrapz(rho,z,Params,Np)

%load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
%load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))
R=Params.R; Rp=Params.Rp;

nodes=linspace(0.00001,0.99999,Np);

Integrand=zeros(length(nodes)); nodesr=Integrand; nodesz=nodesr;
for i=1:length(nodes)

    rhoint=R-Rp+nodes(i)*2*Rp;
    zmax=sqrt(Rp^2-(rhoint-R)^2);
    zmin=-zmax;
    
    for j=1:length(nodes)

        zint=zmin+(zmax-zmin)*nodes(j);
      
        xi=(rho^2+rhoint^2+(z-zint)^2)/(2*rho*rhoint);
        Jphi=Params.J(rhoint,zint);
         
        nodesz(i,j)=zint;
        nodesr(i,j)=rhoint;
        Integrand(i,j)=2*Rp*(zmax-zmin)*sqrt(rhoint/rho)*Q12(xi)*Jphi;

        if Params.test==1
            scatter(rhoint,zint,15)
            hold on
        end
       

    end
end

A=Params.mu/(2*pi)*trapz(nodes,trapz(nodes,Integrand));


end


function A=VectorPotDirMRW(rho,z,Params)

load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))
R=Params.R; Rp=Params.Rp;

Integrand=0;
for i=1:length(nodes)

    rhoip=R+Rp*sqrt(1-(2*nodes(i)-1)^2);
    rhoim=R-Rp*sqrt(1-(2*nodes(i)-1)^2);

    for j=1:length(nodes)
      
        Integrand=Integrand+2*Rp*weights(i)*weights(j)*(rhoip-rhoim)*IntegrandG0(rho,z,rhoim+nodes(j)*(rhoip-rhoim),...
            Rp*(2*nodes(i)-1),Params);

        if Params.test==1
            scatter(rhoint,zint,15)
            hold on
        end
       

    end
end

A=Integrand;


end

function [Bx,Bz]=FieldBDirMRW(rho,z,Params)

load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))
R=Params.R; Rp=Params.Rp;

Integrandx=0; Integrandz=0;
for i=1:length(nodes)

    rhoip=R+Rp*sqrt(1-(2*nodes(i)-1)^2);
    rhoim=R-Rp*sqrt(1-(2*nodes(i)-1)^2);

    for j=1:length(nodes)
      
        [funx,funz]=IntegrandG(rho,z,rhoim+nodes(j)*(rhoip-rhoim),Rp*(2*nodes(i)-1),Params);
        Integrandx=Integrandx+2*Rp*weights(i)*weights(j)*(rhoip-rhoim)*funx;
        Integrandz=Integrandz+2*Rp*weights(i)*weights(j)*(rhoip-rhoim)*funz;

        if Params.test==1
            scatter(rhoint,zint,15)
            hold on
        end
       

    end
end

Bx=Integrandx; Bz=Integrandz;


end

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

function J= Jkid(rho,z,a,b,sigma,R,Rp)
theta=wrapTo2Pi(atan2(z/b,(rho-R)/a));

rad=(rho-R).^2/a^2+z.^2/b^2;

factrad=Rp*(1-0.5*sin(theta/2).^30);

J=exp(-rad.^2./(sigma^2*factrad.^2));

end
function J= Jkid2sign(rho,z,a,b,sigma,R,Rp)
theta=wrapTo2Pi(atan2(z/b,(rho-R)/a));

rad=(rho-R).^2/a^2+z.^2/b^2;

factrad=Rp*(1-0.5*sin(theta/2).^30);

J=exp(-rad.^2./(sigma^2*factrad.^2)).*sin(theta).*rad;

end
function [f,A,B]=PlotfwB(fun,fun2,rhovec,zvec,Params)

A=zeros(length(rhovec),length(zvec)); Bx=A; Bz=A;

vec1=1:length(rhovec);
vec2=1:length(zvec);
parfor i=vec1
    for j=vec2
        
    A(j,i)=fun(rhovec(i),zvec(j));
    [Bx(j,i),Bz(j,i)]=fun2(rhovec(i),zvec(j));
    
    end
end

B=struct('Bx',Bx,'Bz',Bz);
%f=surfc(rhovec,zvec,A,A); hold on
%f2=quiver(rhovec,zvec,Bx,Bz);
%streamline(rhovec,zvec,Bx,Bz,)

% %view(2);
% %shading interp
% %shg
% R=Params.R; Rp=Params.Rp;
% 
% t=linspace(0,2*pi);
% hold on
% plot(R+Rp*cos(t),Rp*sin(t),'r');
% axis equal
% 
% close all

end
function [f,A]=Plotf(fun,rhovec,zvec,Params)

A=zeros(length(rhovec),length(zvec));
for i=1:length(rhovec)
    for j=1:length(zvec)
        
    A(j,i)=fun(rhovec(i),zvec(j));

    end
end

f=surfc(rhovec,zvec,A);

view(2);
shading interp
shg
R=Params.R; Rp=Params.Rp;

t=linspace(0,2*pi);
hold on
plot(R+Rp*cos(t),Rp*sin(t),'r');
axis equal

end

function[Bx,Bz]= BfromA(rho,z,Params)
%Obtaining B from differentiation of A.

delta=Params.DeltaDif;

Bx=-(VectorPot(rho,z+delta,Params)-VectorPot(rho,z-delta,Params))/(2*delta);
Bz=(1/rho)*((rho+delta)*VectorPot(rho+delta,z,Params)-(rho-delta)*...
    VectorPot(rho-delta,z,Params))/(2*delta);
end

function [A]=VectorPot(rho,z,Params)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp; Htheta=Params.Htheta;

if Params.test
    subplot(2,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

alphapunto=atan2(z,rho-R);

if ((rho-R)^2+z^2)<Htheta(alphapunto)^2
    %Psimin=0; Psimax=2*pi;
    Int=zeros(size(nodes)); fun=Int; Funang=fun; Intsum=0;
    for i=1:length(nodes)
        if Params.ExactLims==1
         [~,Dmax]=LimsDIntExact(rho,z,nodes(i)*2*pi,Params);
        else
         [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);
        end

        for j=1:length(nodes)
            fun(j)=IntegrandH(rho,z,nodes(j)*Dmax,nodes(i)*2*pi,Params);
            Int(j)=2*pi*Dmax*fun(j)*weights(i)*weights(j);
        
        end  
        Int(isnan(Int)|isinf(Int))=0;
        Funang(i)=sum(Int);
        Intsum=Intsum+Funang(i);
        if Params.test
            subplot(2,2,2)
            plot(nodes,fun); 
            shg
            %pause();
        end
    end

    if Params.test
            subplot(2,2,4)
            plot(nodes,funang)
            shg
            pause();

    end

else
    if Params.ExactLims==1
        [Psimin,Psimax,Psic,Dc]=LimsPsiExteriorExact(rho,z,Params);
    else
        [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    end
    Int=zeros(size(nodes)); Intsum=0; fun=Int;
    for i=1:length(nodes)
        Psi=Psimin+nodes(i)*(Psimax-Psimin);
        
        if Params.ExactLims==1
            [Dmin,Dmax]=LimsDExtExact(rho,z,Psi,Psimin,Psimax,Params);
        else
            [Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        end

        for j=1:length(nodes)
            D=Dmin+nodes(j)*(Dmax-Dmin);
            fun(j)=IntegrandH(rho,z,D,Psi,Params);
            Int(j)=(Psimax-Psimin)*(Dmax-Dmin)*weights(i)*weights(j)*fun(j);
            %scatter(real(rho+D*cos(Psi)),real(z+D*sin(Psi)),10+10e2*abs(Int(j))); hold on; axis equal
        end       
        Int(isnan(Int)|isinf(Int))=0;
        Intsum=Intsum+sum(Int);
        if Params.test
            subplot(1,2,2)
            plot(nodes,fun)
            pause();
        end
    end


end

A=Intsum;
end

function [A] = VectorPotDE(rho,z,Params,M)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp; 

kvalues=[3,6,12,25,50,101,202,404];
kmax=kvalues(M+1);
kvec=-kmax:kmax;

if Params.test
    subplot(2,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

if ((rho-R)^2+z^2)<=Rp^2
    %Psimin=0; Psimax=2*pi;
    Int=zeros(size(nodes)); fun=Int; Funang=fun; Intsum=0;
    for i=1:length(nodes)
         [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);

        for j=1:length(kvec)
            D=Dmax/2*(1+tanh(pi/2*sinh(kvec(j)/(2^M))));
            fun(j)=IntegrandH(rho,z,D,nodes(i)*2*pi,Params);
            Int(j)=pi^2/(2^(M+1))*Dmax*fun(j)*weights(i)*cosh(kvec(j)/(2^M))/(cosh(pi/2*sinh(kvec(j)/(2^M))))^2;
        
        end  
        Int(isnan(Int)|isinf(Int))=0;
        Funang(i)=sum(Int);
        Intsum=Intsum+Funang(i);
        
        if Params.test
            subplot(2,2,2)
            plot(nodes,fun); 
            shg
            %pause();
        end
    end

    if Params.test
            subplot(2,2,4)
            plot(nodes,funang)
            shg
            pause();

    end

else
    
 %[Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    Int=zeros(size(kvec)); Intsum=0; fun=Int;
    for i=1:length(kvec)
        %Psi=Psimin+nodes(i)*(Psimax-Psimin);
        %[Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        for j=1:length(kvec)
            zj=Rp*tanh(pi/2*sinh(kvec(j)/2^M));
            rhoj2=R+sqrt(Rp^2-zj^2);
            rhoj1=R-sqrt(Rp^2-zj^2);
            %D=Dmin+nodes(j)*(Dmax-Dmin);
            fun(j)=IntegrandG0(rho,z,0.5*(rhoj2-rhoj1)*tanh(pi/2*sinh(kvec(i)/2^M))+0.5*(rhoj2+rhoj1),zj,Params);
            Int(j)=pi^2*Rp/(8*2^(2*M))*(rhoj2-rhoj1)*fun(j)*cosh(kvec(i)/2^M)/(cosh(pi/2*sinh(kvec(i)/2^M))^2)*cosh(kvec(j)/2^M)/(cosh(pi/2*sinh(kvec(j)/2^M))^2);
            %scatter(0.5*(rhoj2-rhoj1)*tanh(pi/2*sinh(kvec(i)/2^M))+0.5*(rhoj2+rhoj1),zj,10+10e2*abs(Int(j))); hold on; %axis equal
        end       
        Int(isnan(Int)|isinf(Int))=0;
        Intsum=Intsum+sum(Int);
        if Params.test
            subplot(1,2,2)
            plot(nodes,fun)
            pause();
        end
    end


end

A=Intsum;


end
function [A]=VectorPotDE6(rho,z,Params,M)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp;
kvalues=[3,6,12,25,50,101,202,404];
kmax=kvalues(M+1);
kvec=-kmax:kmax;


if Params.test
    subplot(2,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

if ((rho-R)^2+z^2)<=Rp^2
    warning('NO USAR VectorPotDE6 En el interior del plasma. Usa E2')
    %Psimin=0; Psimax=2*pi;
    Int=zeros(size(nodes)); fun=Int; Funang=fun; Intsum=0;
    for i=1:length(nodes)
         [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);

        for j=1:length(nodes)
            fun(j)=IntegrandH(rho,z,nodes(j)*Dmax,nodes(i)*2*pi,Params);
            Int(j)=2*pi*Dmax*fun(j)*weights(i)*weights(j);
        
        end  
        Int(isnan(Int)|isinf(Int))=0;
        Funang(i)=sum(Int);
        Intsum=Intsum+Funang(i);
        if Params.test
            subplot(2,2,2)
            plot(nodes,fun); 
            shg
            %pause();
        end
    end

    if Params.test
            subplot(2,2,4)
            plot(nodes,funang)
            shg
            pause();

    end

else
 [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    Int=zeros(size(kvec)); Intsum=0; fun=Int;
    for i=1:length(kvec)
        alphai=(Psimax-Psimin)/2*tanh(pi/2*sinh(kvec(i)/2^M))+(Psimax+Psimin)/2; 
        [Dmin,Dmax]=LimsDExt(alphai,Dc,Psic,Params);
        betaip=Dmax; betaim=Dmin;
        for j=1:length(kvec)
            D=(betaip-betaim)/2*tanh(pi/2*sinh(kvec(j)/2^M))+(betaip+betaim)/2;
            fun(j)=(betaip-betaim)*IntegrandH(rho,z,D,alphai,Params);
            Int(j)=(Psimax-Psimin)*pi^2/(16*2^(2*M))*fun(j)*cosh(kvec(i)/2^M)/(cosh(pi/2*sinh(kvec(i)/2^M))^2)*cosh(kvec(j)/2^M)/(cosh(pi/2*sinh(kvec(j)/2^M))^2);
            %scatter(real(rho+D*cos(alphai)),real(z+D*sin(alphai)),10+10e2*abs(Int(j))); hold on; axis equal
        end       
        Int(isnan(Int)|isinf(Int))=0;
        Intsum=Intsum+sum(Int);
        if Params.test
            subplot(1,2,2)
            plot(nodes,fun)
            pause();
        end
    end


end

A=Intsum;
end
function [A] = VectorPotDE2(rho,z,Params,M)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp; 

kvalues=[3,6,12,25,50,101,202,404];
kmax=kvalues(M+1);
kvec=-kmax:kmax;

if Params.test
    subplot(2,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end


    
 %[Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    Int=zeros(size(nodes)); Intsum=0; fun=Int;
    for i=1:length(kvec)
        %Psi=Psimin+nodes(i)*(Psimax-Psimin);
        %[Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        for j=1:length(kvec)
            zj=Rp*tanh(pi/2*sinh(kvec(j)/2^M));
            rhoj2=R+sqrt(Rp^2-zj^2);
            rhoj1=R-sqrt(Rp^2-zj^2);
            %D=Dmin+nodes(j)*(Dmax-Dmin);
            fun(j)=IntegrandG0(rho,z,0.5*(rhoj2-rhoj1)*tanh(pi/2*sinh(kvec(i)/2^M))+0.5*(rhoj2+rhoj1),zj,Params);
            Int(j)=pi^2*Rp/(8*2^(2*M))*(rhoj2-rhoj1)*fun(j)*cosh(kvec(i)/2^M)/(cosh(pi/2*sinh(kvec(i)/2^M))^2)*cosh(kvec(j)/2^M)/(cosh(pi/2*sinh(kvec(j)/2^M))^2);
            
        end       
        Int(isnan(Int)|isinf(Int))=0;
        Intsum=Intsum+sum(Int);
        if Params.test
            subplot(1,2,2)
            plot(nodes,fun)
            pause();
        end
    end




A=Intsum;


end

function [Bx, Bz]=FieldB(rho,z,Params)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp; Htheta=Params.Htheta;

if Params.test
    subplot(1,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

alphapunto=atan2(z,rho-R);

if ((rho-R)^2+z^2)<=Htheta(alphapunto)^2
    %Psimin=0; Psimax=2*pi;
    Intx=zeros(size(nodes)); funx=Intx; Intsumx=0;
    Intz=zeros(size(nodes)); funz=Intz; Intsumz=0;
    for i=1:length(nodes)
        if Params.ExactLims==1
            [~,Dmax]=LimsDIntExact(rho,z,nodes(i)*2*pi,Params);
        else
            [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);
        end

        for j=1:length(nodes)
            [funx(j),funz(j)]=IntegrandB(rho,z,nodes(j)*Dmax,nodes(i)*2*pi,Params);
            Intx(j)=2*pi*Dmax*funx(j)*weights(i)*weights(j);
            Intz(j)=2*pi*Dmax*funz(j)*weights(i)*weights(j);
        end  
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            shg
            pause();
        end
    end

else
    if Params.ExactLims==1
        [Psimin,Psimax,Psic,Dc]=LimsPsiExteriorExact(rho,z,Params);
    else
        [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    end
    Intx=zeros(size(nodes)); Intsumx=0; funx=Intx;
    Intz=zeros(size(nodes)); Intsumz=0; funz=Intz;
    for i=1:length(nodes)
        Psi=Psimin+nodes(i)*(Psimax-Psimin);
        if Params.ExactLims==1
            [Dmin,Dmax]=LimsDExtExact(rho,z,Psi,Psimin,Psimax,Params);
        else
            [Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        end
        for j=1:length(nodes)
            D=Dmin+nodes(j)*(Dmax-Dmin);
            [funx(j),funz(j)]=IntegrandB(rho,z,D,Psi,Params);
            Intx(j)=(Psimax-Psimin)*(Dmax-Dmin)*weights(i)*weights(j)*funx(j);
            Intz(j)=(Psimax-Psimin)*(Dmax-Dmin)*weights(i)*weights(j)*funz(j);

        end       
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            pause();
        end
    end


end
Bx=Intsumx; Bz=Intsumz;

%B=struct('Bx',Bx,'Bz',Bz);

end
function [Bx, Bz]=FieldBDE(rho,z,Params,M)
load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW))
load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW))

R=Params.R; Rp=Params.Rp;
kvalues=[3,6,12,25,50,101,202,404];
kmax=kvalues(M+1);
kvec=-kmax:kmax;

if Params.test
    subplot(1,2,1)
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
hold on; axis equal; plot(rho,z,'o')
%pause()
end

if ((rho-R)^2+z^2)<=Rp^2
    %Psimin=0; Psimax=2*pi;
    Intx=zeros(size(nodes)); funx=Intx; Intsumx=0;
    Intz=zeros(size(nodes)); funz=Intz; Intsumz=0;
    for i=1:length(nodes)
         [~,Dmax]=LimsDInt(rho,z,nodes(i)*2*pi,Params);

        for j=1:length(kvec)
            D=Dmax/2*(1+tanh(pi/2*sinh(kvec(j)/(2^M))));
            [funx(j),funz(j)]=IntegrandB(rho,z,D,nodes(i)*2*pi,Params);
            Intx(j)=pi^2/(2^(M+1))*Dmax*funx(j)*weights(i)*cosh(kvec(j)/(2^M))/(cosh(pi/2*sinh(kvec(j)/(2^M))))^2;
            Intz(j)=pi^2/(2^(M+1))*Dmax*funz(j)*weights(i)*cosh(kvec(j)/(2^M))/(cosh(pi/2*sinh(kvec(j)/(2^M))))^2;
        end  
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            shg
            pause();
        end
    end

else
 %[Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params);
    Intx=zeros(size(nodes)); Intsumx=0; funx=Intx;
    Intz=zeros(size(nodes)); Intsumz=0; funz=Intz;
    for i=1:length(kvec)
        %Psi=Psimin+nodes(i)*(Psimax-Psimin);
        %[Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params);
        for j=1:length(kvec)
            zj=Rp*tanh(pi/2*sinh(kvec(j)/2^M));
            rhoj2=R+sqrt(Rp^2-zj^2);
            rhoj1=R-sqrt(Rp^2-zj^2);
            [funx(j),funz(j)]=IntegrandG(rho,z,0.5*(rhoj2-rhoj1)*tanh(pi/2*sinh(kvec(i)/2^M))+0.5*(rhoj2+rhoj1),zj,Params);
            Intx(j)=pi^2*Rp/(8*2^(2*M))*(rhoj2-rhoj1)*funx(j)*cosh(kvec(i)/2^M)/(cosh(pi/2*sinh(kvec(i)/2^M))^2)*cosh(kvec(j)/2^M)/(cosh(pi/2*sinh(kvec(j)/2^M))^2);
            Intz(j)=pi^2*Rp/(8*2^(2*M))*(rhoj2-rhoj1)*funz(j)*cosh(kvec(i)/2^M)/(cosh(pi/2*sinh(kvec(i)/2^M))^2)*cosh(kvec(j)/2^M)/(cosh(pi/2*sinh(kvec(j)/2^M))^2);

        end       
        Intx(isnan(Intx)|isinf(Intx))=0;
        Intsumx=Intsumx+sum(Intx);
        Intz(isnan(Intz)|isinf(Intz))=0;
        Intsumz=Intsumz+sum(Intz);
        if Params.test
            subplot(1,2,2)
            plot(nodes,funx)
            pause();
        end
    end


end
Bx=Intsumx; Bz=Intsumz;

%B=struct('Bx',Bx,'Bz',Bz);

end
function [Dmin,Dmax]=LimsDExt(Psi,Dc,Psic,Params)
Rp=Params.Rp;
Dmax=Dc*abs(cos(Psi-Psic))+sqrt(Rp^2-Dc^2*sin(Psi-Psic)^2);
Dmin=Dc*abs(cos(Psi-Psic))-sqrt(Rp^2-Dc^2*sin(Psi-Psic)^2);

end

function [Psimin,Psimax,Psic,Dc]=LimsPsiExterior(rho,z,Params)
R=Params.R; Rp=Params.Rp;
Psic=atan2(-z,(R-rho));
Dc=sqrt((R-rho)^2+z^2);

Psimax=Psic+asin(Rp/Dc);
Psimin=Psic-asin(Rp/Dc);

end


function [Dmin,Dmax]=LimsDInt(rho,z,Psi,Params)
R=Params.R; Rp=Params.Rp;

Dmin=0;
Dmax=((R-rho)*cos(Psi)-z*sin(Psi))+sqrt(((R-rho)*cos(Psi)-z*sin(Psi))^2+(Rp^2-z^2-(R-rho)^2));

end

function [Dmin,Dmax]=LimsDIntExact(rho,z,Psi,Params)
R=Params.R; Rp=Params.Rp; Htheta=Params.Htheta;


Dmin=0;
minfun=@(alpha) abs(cos(Psi)*(Htheta(alpha)*sin(alpha)-z)-sin(Psi)*(Htheta(alpha)*cos(alpha)+(R-rho)));
alphap=Psi+asin((sin(Psi)*(R-rho)+z*cos(Psi))/Rp);
alphasol=fminsearch(minfun,alphap);
Dmax=((R-rho)*cos(Psi)-z*sin(Psi))+sqrt(((R-rho)*cos(Psi)-z*sin(Psi))^2+(Htheta(alphasol)^2-z^2-(R-rho)^2));

end

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
while abs(wrapToPi(alphasol2-alphasol1))<1e-1 || feval2>0.01
[alphasol2,feval2]=fminsearch(minfun,alphaprop2);
alphaprop2=alphaprop2+rand*pi; n2=n2+1;
if abs(wrapToPi(Psimin-Psi))<0.3 || abs(wrapToPi(Psimax-Psi))<0.3
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

function [H1,H2]=IntegrandB(rho,z,D,Psi,Params)


rhop=rho+D*cos(Psi); zp=z+D*sin(Psi);
if Params.test
    subplot(1,2,1)
plot(rhop,zp,'x'); shg
%pause();
end

xi=(rho^2+rhop^2+(z-zp)^2)/(2*rho*rhop);

DerQhalf=DerQ12(xi);
Qhalf=Q12(xi);

Jphi=Params.J(rhop,zp);

H1=-Params.mu/(2*pi)*(z-zp)/(sqrt(rhop)*rho^(3/2))*DerQhalf*Jphi*D;

H2=Params.mu/(2*pi)*(sqrt(rhop)/(2*rho^(3/2))*Qhalf+...
    ((rho^2-rhop^2-(z-zp)^2)/(2*rho^(5/2)*sqrt(rhop)))*DerQhalf)*Jphi*D;

end


function H0=IntegrandH(rho,z,D,Psi,Params)


rhop=rho+D*cos(Psi); zp=z+D*sin(Psi);
if Params.test
    subplot(2,2,1)
    hold on;
plot(rhop,zp,'x'); shg
%pause();
end

xi=(rho^2+rhop^2+(z-zp)^2)/(2*rho*rhop);

[Qhalf,errorc]=Q12(xi);

if errorc&&Params.test
plot(rhop,zp,'x'); hold on;
plot(rho,z,'o');
R=Params.R; Rp=Params.Rp;
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
axis equal
shg
pause();
end

Jphi=Params.J(rhop,zp);

H0=Params.mu/(2*pi)*sqrt(rhop/rho)*Qhalf*Jphi*D;

end

function [A, Bx, Bz]=FieldLoop(rho,z,Params)
a=Params.R; mu0=Params.mu;
I=1;

%A=mu0*I*sqrt(R)/(2*pi*sqrt(rho))*Q12((rho^2+a^2+z^2)/2*a*rho);

k2=4*a*rho/((rho+a)^2+z^2);
[K,E]=ellipke(k2);

A=mu0*I*sqrt(a)/(2*pi*sqrt(rho))*((2-k2)*K-2*E);

Bx=mu0*I*z/(2*pi*rho*sqrt((rho+a)^2+z^2))*(-K+(rho^2+a^2+z^2)/((rho-a)^2+z^2)*E);

Bz=mu0*I/(2*pi*sqrt((rho+a)^2+z^2))*(K+(a^2-rho^2-z^2)/((rho-a)^2+z^2)*E);

end

function [A, Bx, Bz]= FieldCilinder(rho,z,Params)

theta=atan2(z,rho-Params.R);
rad=sqrt(z^2+(rho-Params.R)^2);

if rad<Params.Rp
Bo=Params.mu*rad/2;
Bx=Bo*sin(theta); Bz=-Bo*cos(theta);
A=-Params.mu*rad^2/4;
else
Bo=Params.mu*Params.Rp^2/(2*rad);
Bx=Bo*sin(theta); Bz=-Bo*cos(theta);
A=0.5*Params.mu*Params.Rp^2*(log(Params.Rp/rad)-0.5);
end




end

function DerQhalf= DerQ12(x)
% if imag(x)~=0
%     a=1;
% end

%[K,E]=ellipke(sqrt(2./(x+1)));
try
[K,E]=ellipke(2/(x+1));
DerQhalf=0.5*sqrt(2/(x+1))*K-x*(2*(x+1))^0.5/(2*(x^2-1))*E;
catch
DerQhalf=0; errorc=1;
end
end
function [Qhalf,errorc]= Q12(x)


% if imag(x)~=0
%     a=1;
% end
M=(2/(x+1));



try
%[K,E]=ellipke(sqrt(2./(x+1)));
[K,E]=ellipke(M);
Qhalf=x*sqrt(2/(x+1))*K-sqrt(2*(x+1))*E;
errorc=0;
catch
    %warning('error en q1/2')
    errorc=1; Qhalf=0;
end

end


function G0=IntegrandG0(rho,z,rhop,zp,Params)
% Integrand of A in cylindrical coordinates

%rhop=rho+D*cos(Psi); zp=z+D*sin(Psi);
if Params.test
    subplot(2,2,1)
    hold on;
plot(rhop,zp,'x'); shg
%pause();
end

xi=(rho^2+rhop^2+(z-zp)^2)/(2*rho*rhop);

[Qhalf,errorc]=Q12(xi);

if errorc&&Params.test
plot(rhop,zp,'x'); hold on;
plot(rho,z,'o');
R=Params.R; Rp=Params.Rp;
plot(R+Rp*cos(linspace(0,2*pi)),Rp*sin(linspace(0,2*pi)));
axis equal
shg
pause();
end

Jphi=Params.J(rhop,zp);

G0=Params.mu/(2*pi)*sqrt(rhop/rho)*Qhalf*Jphi;

end

function [G1,G2]=IntegrandG(rho,z,rhop,zp,Params)
% Integrands of B in cylindrical coordinates

%rhop=rho+D*cos(Psi); zp=z+D*sin(Psi);
if Params.test
    subplot(1,2,1)
plot(rhop,zp,'x'); shg
%pause();
end

xi=(rho^2+rhop^2+(z-zp)^2)/(2*rho*rhop);

DerQhalf=DerQ12(xi);
Qhalf=Q12(xi);

Jphi=Params.J(rhop,zp);

G1=-Params.mu/(2*pi)*(z-zp)/(sqrt(rhop)*rho^(3/2))*DerQhalf*Jphi;

G2=Params.mu/(2*pi)*(sqrt(rhop)/(2*rho^(3/2))*Qhalf+...
    ((rho^2-rhop^2-(z-zp)^2)/(2*rho^(5/2)*sqrt(rhop)))*DerQhalf)*Jphi;

end
