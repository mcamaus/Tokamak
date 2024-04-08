function [ B ] = eval_B( r )
%Devuelve las componentes de B considerando simetría en la dirección
%poloidal

global mu0 Xmax Xmin Zmax Zmin Nx Nz Bx By Bz
% global  Nz Nr Rmax Rmin Zmax Zmin BPolaresVertical BPolaresToroidal BPolaresRadial

% Se obtienen las coordenadas polares y relativas del punto
if r(1)>0
    PitchParticula = atan (r(2)/r(1)) ;
else
    PitchParticula = atan (r(2)/r(1)) + pi ;
end

if PitchParticula < 0
PitchParticula = PitchParticula+2*pi;
end

RadioParticula = sqrt(r(1)^2 + r(2)^2);

% Se interpola en la malla inicial y se obtienen las componentes del campo
% en r
Bradial = interp2(linspace(Xmin,Xmax,Nx),linspace(Zmin,Zmax,Nz),Bx,RadioParticula,r(3));
Btoroidal  = interp2(linspace(Xmin,Xmax,Nx),linspace(Zmin,Zmax,Nz),By,RadioParticula,r(3));
Bzeta = interp2(linspace(Xmin,Xmax,Nx),linspace(Zmin,Zmax,Nz),Bz,RadioParticula,r(3));

% Se calcula B

B(1,1) = Bradial*cos(PitchParticula)-Btoroidal*sin(PitchParticula); 
B(2,1) = Bradial*sin(PitchParticula)+Btoroidal*cos(PitchParticula); 
B(3,1) = Bzeta; 

end

