%% Fitting of D parameters to nt_data.mat

%% Fig 9a

R=0.5; Rp=0.42;


%R=0.5; Rp=0.001;

%Currentshape="Constant";
%Currentshape="Gaussian"; 
%Currentshape="Kidney"; 
%Currentshape="2SignKid"; 
Currentshape="invD";


load("nt_data.mat")
j=j*1e-6;

Jrhoz=@(rhop,zp) interp2(r,z,j,rhop,zp);

Params=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',15,'DeltaDif',1e-6, 'test',0);

Params.MRWn=load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params.MRW)).nodes;
Params.MRWw=load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params.MRW)).weights;

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

%% 
%figure
[M,c]=contour(xmat,zmat,Jrhom,linspace(1e-2,1e-1,5),"ShowText",true);


level=M(1,1); Npoints=M(2,1);
xlevel=M(1,2:1+Npoints);
ylevel=M(2,2:1+Npoints);

%hold on
plot(xlevel,ylevel);


%%
R=0.465; Rp=0.4;
delta=-0.3; beta=asin(delta); k=1.7;

rhof=@(alpha) R+Rp/k*cos(alpha+beta*sin(alpha));
zf=@(alpha) Rp*sin(alpha);
alpha=linspace(-pi,pi);
dists=sqrt((rhof(alpha)-R).^2+zf(alpha).^2);
maxdist=max(dists);
Rp=maxdist;


xvec=linspace(0.1,R+5*0.15,100); x=xvec;
zvec=linspace(-5*0.15,5*0.15,100); z=zvec;


%figure;
hold on
plot3(rhof(alpha),zf(alpha),12*ones(size(alpha)),'-','LineWidth',2)
