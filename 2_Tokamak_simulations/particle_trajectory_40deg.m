 close all
clc
% Primero hay que ejecutar GeneraB_Dcoil.m

global q m D a Emax

% Integration final time
t = 10E-5;
D=.25;
a=.125;
q = 1.6E-19; %Particle charge (proton)
m = 1.67E-27*3; %Particle mass (proton)
% Particle initial energy (eV)
part_energ_eV=1e3; 
part_energ_0 = part_energ_eV*1.602E-19; 

% % % Valor máximo del campo eléctrico radial (En V/m)
% Emax = -380E3;

% Particle Initial conditions
r0 = [1.25*R 0 0]';
B_r0=eval_B(r0);
B_ro_unit=B_r0/norm(B_r0);
vect_perp_2_B=cross(B_ro_unit,[1 0 0]);
vect_perp_2_B_unit=vect_perp_2_B/norm(vect_perp_2_B);
pitch_angle=40*pi/180;  % I impose a pitch angle
v_parallel=cos(pitch_angle)*B_ro_unit;
v_perp=sin(pitch_angle)*vect_perp_2_B_unit';
v0_unit=v_parallel+v_perp;
v0 = v0_unit * sqrt( part_energ_0 * 2 / m );

pitch_angle_deg=acos(B_ro_unit'*v0_unit)*180/pi

% Integrador numérico para la particula
tic
options = odeset('Reltol',1e-5);
[T, Trayectoria] = ode45 (@lorentz_law, [0 t], [r0 ; v0],options);
Trayectoria = Trayectoria';
time_to_integrate = toc
plot_traject_tokamak_MCA


