% Integrador num?rico para seguir lineas de campo
tic
close all
clc

% Primero hay que ejecutar GeneraBgrid_2.m y plot_B_tokamak
%
% Integration final time
clear traject
% list_field_line=[];
t_int = 200;
% Mesh for poincare plot
%mesh scale wrt plasma radii
poincare_mesh_scale=1;
n_lines_x=3;
x_min=plasma_center_xcoord-r_plasma*poincare_mesh_scale;
x_max=plasma_center_xcoord+r_plasma*poincare_mesh_scale;
n_lines_z=3;
z_min=-r_plasma*poincare_mesh_scale;
z_max=r_plasma*poincare_mesh_scale;
[x_mesh_poincare,z_mesh_poincare]=meshgrid(linspace(x_min,x_max,n_lines_x),linspace(z_min,z_max,n_lines_z));
y_mesh_poincare=zeros(size(x_mesh_poincare));

for i=1:n_lines_z
    for j=1:n_lines_x
        r0_mesh_point = [x_mesh_poincare(i,j) y_mesh_poincare(i,j) z_mesh_poincare(i,j)]';
        options = odeset('Reltol',1e-4);
        [T, field_line] = ode45 (@field_line_eq, [0 t_int], r0_mesh_point,options);
        traject(i,j).time=T;
        traject(i,j).coords=field_line;
        i,j
    end
end
toc
%}

% Esto es para pintar el poincare plot
field_line_tracing_plot












