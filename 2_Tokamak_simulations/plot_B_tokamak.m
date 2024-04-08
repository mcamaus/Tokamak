% Figuras para analizar el campo magn?tico, run generaB_Dcoils before
global plasma_center_xcoord r_plasma
figure(1)
hold on
plot3([0 0],[0 0],[-.2 .2],'b')
for i=1:N_toroidal_coils
plot3(toroidal_coil_data(i).coords_coil(:,1),toroidal_coil_data(i).coords_coil(:,2)...
    ,toroidal_coil_data(i).coords_coil(:,3),'r')
end
plot3(toroidal_coil_data(1).coords_coil(:,1),toroidal_coil_data(1).coords_coil(:,2)...
    ,toroidal_coil_data(1).coords_coil(:,3),'b')
% plot3(plasma_current_data(1).coords_coil(:,1),plasma_current_data(1).coords_coil(:,2)...
%     ,plasma_current_data(1).coords_coil(:,3),'b')
plot3(lower_BvL_coil_data(1).coords_coil(:,1),lower_BvL_coil_data(1).coords_coil(:,2)...
    ,lower_BvL_coil_data(1).coords_coil(:,3),'K')
plot3(upper_BvL_coil_data(1).coords_coil(:,1),upper_BvL_coil_data(1).coords_coil(:,2)...
    ,upper_BvL_coil_data(1).coords_coil(:,3),'K')
plot3(lower_BvU_coil_data(1).coords_coil(:,1),lower_BvU_coil_data(1).coords_coil(:,2)...
    ,lower_BvU_coil_data(1).coords_coil(:,3),'K')
plot3(upper_BvU_coil_data(1).coords_coil(:,1),upper_BvU_coil_data(1).coords_coil(:,2)...
    ,upper_BvU_coil_data(1).coords_coil(:,3),'K')

% plot3(r_evaluate_B(1),r_evaluate_B(2),r_evaluate_B(3),'xk')
xlabel('X coordenate (m)')
ylabel('Y coordenate (m)')
zlabel('Z coordenate (m)')

view(35,30)
axis([-1,1,-1,1,-1,1])
box
grid
hold off

% figure(2)
% hold on
% plot(1:N_toroidal_coils,[toroidal_coil_data.mod_B],'b')
% xlabel('coil index')
% ylabel('magnitude of B (T)')
% 
% % view(35,30)
% % axis([-3,3,-3,3,-3,3])
% % box
% % grid
% % hold off
% B_max=max(max(mod_B_mesh));
% 

% figure(3)
figure('pos',[600 100 460 400*1.1/.55])
hold on
set(gca,'Color','k');
omega_II_colormap
% "phong" lighting is good for curved, interpolated surfaces. "gouraud"
% is also good for curved surfaces
surf(x_mesh,z_mesh,mod_B_mesh,'AmbientStrength',1);

% quiver(Bx,Bz,10)
view(0,90);
% view(35,40);
plot3(toroidal_coil_data(1).coords_coil(:,1),toroidal_coil_data(1).coords_coil(:,3)...
    ,100.*10*ones(length(toroidal_coil_data(1).coords_coil(:,2))),'k')
plot3(plasma_center_xcoord,0,10,'xk','LineWidth',1.5)
plot3(R_upper_BvL_coil,h_upper_BvL_coil,10,'xk','LineWidth',1.5)
plot3(R_lower_BvL_coil,h_lower_BvL_coil,10,'xk','LineWidth',1.5)
plot3(R_upper_BvU_coil,h_upper_BvU_coil,10,'xk','LineWidth',1.5)
plot3(R_lower_BvU_coil,h_lower_BvU_coil,10,'xk','LineWidth',1.5)
plot3(R_upp_DIV_coil,h_upp_DIV_coil,10,'xk','LineWidth',1.5)
plot3(R_low_DIV_coil,h_low_DIV_coil,10,'xk','LineWidth',1.5)

plot3(vessel_profile_coords(1,:),vessel_profile_coords(2,:),10*ones(length(vessel_profile_coords)),'-k','LineWidth',1.5)
%plot3(plasma_profile(:,1),plasma_profile(:,2),10*ones(length(plasma_profile)),'--r','LineWidth',1.5)

axis([0,.6,-.7,.7])
xlabel('radial coordinate (m)');
ylabel('vertical coordinate (m)');
c = colorbar;
ylabel(c, 'Field modulus (T)');
caxis([0,.8])
shading interp;
lighting gouraud;

%%
figure

streamslice(x_mesh,z_mesh,Bxp,Bzp);




