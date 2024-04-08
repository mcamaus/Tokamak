 %% Tokamak: D-shaped toroidal coils, 2 poloidal coils and plasma current field.
tic
clear all, clc;
global mu0 Xmax Xmin Zmax Zmin Nx Nz Bx By Bz plasma_center_xcoord r_plasma
%% Parameters definition
%% Tokamak parameters
R_tok=.25; % Major radius of the torus
mu0=4*pi*1E-7;% Vacuum permeability
plasma_center_xcoord=0.25;
%vessel local coordinates
vessel_profile_x_local_coord=[0 .26 .36 .36 .26 0 0];
vessel_profile_y_local_coord=[.32 .32 .16 -.16 -.32 -.32 .32];
%displacement to put vessel in place
vessel_profile_x_coord=vessel_profile_x_local_coord+.104;
vessel_profile_y_coord=vessel_profile_y_local_coord;
vessel_profile_coords=[vessel_profile_x_coord ; vessel_profile_y_coord];


%%                                          Coils
%                                         POLOIDAL FIELD COILS
%             PLASMA: Modeled as a circular coil
N_segments_plasma_current=30;   % Number of segments of plama current:
h_plasma=0;                     % Plasma vertical position
I_plasma=-50e3;                   % Plasma current (A)
% I_plasma=0;                   % Plasma current (A)
center_plasma_current=[0 0 h_plasma];
rotations_plasma_current=[0 0 0];

%             PLASMA: Modeled as a distributed current density
r_plasma=.42/2;
% Plasma boundary for plasma current
alpha=0:.1:2*pi;
plasma_profile=[cos(alpha); sin(alpha)]'*r_plasma;
plasma_profile(:,1)=plasma_profile(:,1)+plasma_center_xcoord;

%               Parameters for plasma's B field calculation  

Rp=r_plasma; R=plasma_center_xcoord; Currentshape="Constant";
%only constant current is currently normalized
R=0.5/2; Rp=0.42/2;


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

%j=j;

Jrhoz=@(rhop,zp) 1.5*interp2(r/2,z/2,j,rhop,zp);



end


Params_Plasma=struct('R',R, 'Rp', Rp,'J',Jrhoz,'mu',4*pi*1e-7,'MRW',15,'DeltaDif',1e-6, 'test',0);
Params_Plasma.MRWn=load(sprintf('MRW_quadrature/NodesMRW%2i.mat',Params_Plasma.MRW)).nodes;
Params_Plasma.MRWw=load(sprintf('MRW_quadrature/WeightsMRW%2i.mat',Params_Plasma.MRW)).weights;

%%              Lower (BvL) UPPER COIL 
N_segments_upp_BvL_coil=20; % Number of segments of upper poloidal coil:
R_upper_BvL_coil=.5475;           % Poloidal upper coil radius (m)
h_upper_BvL_coil=.17;            % Poloidal upper coil vertical position (m)
I_upper_BvL_coil=-5e3;         % Poloidal upper coil current (amperes)
% I_upper_BvL_coil=0;         % Poloidal upper coil current (amperes)
center_upper_BvL_coil=[0 0 h_upper_BvL_coil];
rotations_upper_BvL_coil=[0 0 0];
%%              Lower (BvL) LOWER COIL
N_segments_low_BvL_coil=20; % Number of segments of lower poloidal coil:
R_lower_BvL_coil= .5475;          % Poloidal lower coil radius (m)
h_lower_BvL_coil=-.17;           % Poloidal lower coil vertical position (m)
I_lower_BvL_coil=-5e3;         % Poloidal lower coil current (amperes)
% I_lower_BvL_coil=0;         % Poloidal upper coil current (amperes)
center_lower_BvL_coil=[0 0 h_lower_BvL_coil];
rotations_lower_BvL_coil=[0 0 0];
%%              Upper (BvU)  UPPER COIL
N_segments_upp_BvU_coil=20; % Number of segments of upper poloidal coil:
R_upper_BvU_coil=.44;           % Poloidal upper coil radius (m)
h_upper_BvU_coil=.416;            % Poloidal upper coil vertical position (m)
I_upper_BvU_coil=5e3;         % Poloidal upper coil current (amperes)
% I_upper_BvU_coil=0;         % Poloidal upper coil current (amperes)
center_upper_BvU_coil=[0 0 h_upper_BvU_coil];
rotations_upper_BvU_coil=[0 0 0];
%%              Upper (BvU)  LOWER COIL
N_segments_low_BvU_coil=20; % Number of segments of lower poloidal coil:
R_lower_BvU_coil= .44;          % Poloidal lower coil radius (m)
h_lower_BvU_coil=-.416;           % Poloidal lower coil vertical position (m)
I_lower_BvU_coil=5e3;         % Poloidal lower coil current (amperes)
% I_lower_BvU_coil=0;         % Poloidal upper coil current (amperes)
center_lower_BvU_coil=[0 0 h_lower_BvU_coil];
rotations_lower_BvU_coil=[0 0 0];

%%              Upper DIVERTOR COIL (DIV)
N_segments_upp_DIV_coil=20; % Number of segments of lower poloidal coil:
R_upp_DIV_coil= .22;          % Poloidal lower coil radius (m)
h_upp_DIV_coil=.34;           % Poloidal lower coil vertical position (m)
I_upp_DIV_coil=30e3;         % Poloidal lower coil current (amperes)
% I_upp_DIV_coil=0;         % Poloidal upper coil current (amperes)
center_upp_DIV_coil=[0 0 h_upp_DIV_coil];
rotations_upp_DIV_coil=[0 0 0];

%%              Lower DIVERTOR COIL (DIV)
N_segments_low_DIV_coil=20; % Number of segments of lower poloidal coil:
R_low_DIV_coil= .22;          % Poloidal lower coil radius (m)
h_low_DIV_coil=-.34;           % Poloidal lower coil vertical position (m)
I_low_DIV_coil=30e3;         % Poloidal lower coil current (amperes)
% I_low_DIV_coil=0;         % Poloidal upper coil current (amperes)
center_low_DIV_coil=[0 0 h_low_DIV_coil];
rotations_low_DIV_coil=[0 0 0];


%%                                          TOROIDAL FIELD COILS
N_toroidal_coils=8;        % Number of toroidal coils
geom_Dshape_coil;  %Code to generate Dshape coil profile coordinates
r_coil_center= 0.042; %Dist from coil local center to tokamak center
coords_local_coil=[Dcoil_profile_x_coords' Dcoil_profile_y_coords' zeros(1,length(Dcoil_profile_x_coords))'];
I_toroidal_coil=60E3;      % Toroidal coil current (A)
% I_toroidal_coil=0;      % Toroidal coil current (A)
angle_toroidal_coil=0:2*pi/N_toroidal_coils:2*pi;   %Angle of each coil
center_toroidal_coil=[r_coil_center*cos(angle_toroidal_coil)' r_coil_center*sin(angle_toroidal_coil)'...
    0*ones(1,length(angle_toroidal_coil))'];        %Center of each coil

%% Evaluate magnetic field at each point of a mesh in space
% 3D cartesian coordinates of the point where B is to be evaluated
%% 3D mesh to evaluate B at each node
Xmin = 0 ;        %Minimum radius
Xmax = max(vessel_profile_x_coord)*1.3 ; %Maximum radius
Zmin = -max(coords_local_coil(:,2))*1.1 ;%Maximum vertical position
Zmax =  max(coords_local_coil(:,2))*1.1 ;%Minimum vertical position
Nx=157; %Number of points in radial direction for the mesh
Nz=157; %Number of points in vertical direction for the mesh
[x_mesh,z_mesh]=meshgrid(linspace(Xmin,Xmax,Nx),linspace(Zmin,Zmax,Nz));
y_mesh=zeros(size(x_mesh));
lista_NaN=[];

%Bxp=zeros(Nx,Nz);


for k = Nz:-1:1
    % for k = 7
    
    for j = Nx:-1:1
        %     for j = 13
        r_evaluate_B = [x_mesh(j,k) y_mesh(j,k) z_mesh(j,k)]';
        %% Manetic field produced by toroidal coils
        for i=1:N_toroidal_coils
            %         for i=1
            % 3 rotations to orient coil
            rotations=[pi/2 0 2*pi/N_toroidal_coils*(i-1)];
%             [toroidal_coil_data(i)]=magn_field_circ_coil_3D(I_toroidal_coil,r_evaluate_B,...
%                 center_toroidal_coil(i,:),R_toroidal_coil,N_segments_toroidal_coil,rotations);  %Coil circular            
            [toroidal_coil_data(i)]=magn_field_Dshape_coil_3D(I_toroidal_coil,r_evaluate_B,...
                center_toroidal_coil(i,:),coords_local_coil,rotations);  %Coil circular
        end
        % Data for each point where B is evaluated
        % Add the field generated by each coil
        mesh_point_data(j,k).B=sum([toroidal_coil_data.B],2);
        mesh_point_data(j,k).coords=r_evaluate_B;
        mesh_point_data(j,k).toroidal_coil_data=toroidal_coil_data;
        
        %% Magnetic field produced by plasma current
        % Plasma modeled as a circular coil
        plasma_center_radius=plasma_center_xcoord;
        %[plasma_current_data]=magn_field_circ_coil_3D(I_plasma,r_evaluate_B,...
         %   center_plasma_current,plasma_center_radius,N_segments_plasma_current,rotations_plasma_current);
        %if isnan(plasma_current_data.mod_B)
        %   lista_NaN=[lista_NaN; j k]; 
        %end
        % Plasma modeled as a distributed current density in toroid
%         [plasma_current_data]=magn_field_distrib_plasma_toroid_3D(I_plasma,r_evaluate_B,...
%             center_plasma_current, plasma_center_xcoord, r_plasma);
        
        [plasma_current_data]=FieldB_plasma(r_evaluate_B,Params_Plasma);
        
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+plasma_current_data.B;
        mesh_point_data(j,k).Bp=plasma_current_data.B;
        mesh_point_data(j,k).plasma_current_data=plasma_current_data;
        
        %% Magnetic field produced by BvL upper poloidal coil
        [upper_BvL_coil_data]=magn_field_circ_coil_3D(I_upper_BvL_coil,r_evaluate_B,...
            center_upper_BvL_coil,R_upper_BvL_coil,N_segments_upp_BvL_coil,rotations_upper_BvL_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+upper_BvL_coil_data.B;
        mesh_point_data(j,k).upper_BvL_coil_data=upper_BvL_coil_data;
        
        %% Magnetic field produced by BvL lower poloidal coil
        [lower_BvL_coil_data]=magn_field_circ_coil_3D(I_lower_BvL_coil,r_evaluate_B,...
            center_lower_BvL_coil,R_lower_BvL_coil,N_segments_low_BvL_coil,rotations_lower_BvL_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+lower_BvL_coil_data.B;
        mesh_point_data(j,k).lower_BvL_coil_data=lower_BvL_coil_data;
        
        %% Magnetic field produced by BvU upper poloidal coil
        [upper_BvU_coil_data]=magn_field_circ_coil_3D(I_upper_BvU_coil,r_evaluate_B,...
            center_upper_BvU_coil,R_upper_BvU_coil,N_segments_upp_BvU_coil,rotations_upper_BvU_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+upper_BvU_coil_data.B;
        mesh_point_data(j,k).upper_BvU_coil_data=upper_BvU_coil_data;
        
        %% Magnetic field produced by BvU lower poloidal coil
        [lower_BvU_coil_data]=magn_field_circ_coil_3D(I_lower_BvU_coil,r_evaluate_B,...
            center_lower_BvU_coil,R_lower_BvU_coil,N_segments_low_BvU_coil,rotations_lower_BvU_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+lower_BvU_coil_data.B;
        mesh_point_data(j,k).lower_BvU_coil_data=lower_BvU_coil_data;
      
        %% Magnetic field produced by DIV upper poloidal coil
        [upper_DIV_coil_data]=magn_field_circ_coil_3D(I_upp_DIV_coil,r_evaluate_B,...
            center_upp_DIV_coil,R_upp_DIV_coil,N_segments_upp_DIV_coil,rotations_upp_DIV_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+upper_DIV_coil_data.B;
        mesh_point_data(j,k).upper_DIV_coil_data=upper_DIV_coil_data;
        
        %% Magnetic field produced by DIV upper poloidal coil
        [lower_DIV_coil_data]=magn_field_circ_coil_3D(I_low_DIV_coil,r_evaluate_B,...
            center_low_DIV_coil,R_low_DIV_coil,N_segments_low_DIV_coil,rotations_low_DIV_coil);
        mesh_point_data(j,k).B=mesh_point_data(j,k).B+lower_DIV_coil_data.B;
        mesh_point_data(j,k).lower_DIV_coil_data=lower_DIV_coil_data;

        %% Magnetic field magnitude
        mod_B_mesh(j,k)=sqrt(mesh_point_data(j,k).B(1)^2+...
            mesh_point_data(j,k).B(2)^2+mesh_point_data(j,k).B(3)^2);
        
        mod_Bp_mesh(j,k)=sqrt(mesh_point_data(j,k).Bp(1)^2+...
            mesh_point_data(j,k).Bp(2)^2+mesh_point_data(j,k).Bp(3)^2);
        %% Magnetic field components 
        Bx(j,k)=mesh_point_data(j,k).B(1);
%                 Bx(j,k)=zeros(size(mesh_point_data(j,k).B(1)));

        By(j,k)=mesh_point_data(j,k).B(2);
        Bz(j,k)=mesh_point_data(j,k).B(3);

        Bxp(j,k)=mesh_point_data(j,k).Bp(1);
%                 Bx(j,k)=zeros(size(mesh_point_data(j,k).B(1)));

        Byp(j,k)=mesh_point_data(j,k).Bp(2);
        Bzp(j,k)=mesh_point_data(j,k).Bp(3);
        
%                 Bz(j,k)=zeros(size(mesh_point_data(j,k).B(3)));

    end
end
%         %% Interpolate singularities
%         Bx(lista_NaN(1),lista_NaN(2))=(Bx(lista_NaN(1)-1,lista_NaN(2))+Bx(lista_NaN(1)+1,lista_NaN(2)))/2;
%         By(lista_NaN(1),lista_NaN(2))=(By(lista_NaN(1)-1,lista_NaN(2))+By(lista_NaN(1)+1,lista_NaN(2)))/2;
%         Bz(lista_NaN(1),lista_NaN(2))=(Bz(lista_NaN(1)-1,lista_NaN(2))+Bz(lista_NaN(1)+1,lista_NaN(2)))/2;
toc

plot_B_tokamak;
% field_line_tracing;






