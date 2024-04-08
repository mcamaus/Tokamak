% Plot 3D trajectory+vessel
%
% figure(1)

vessel_sector_ini=45*pi/180;
vessel_sector_end=180*pi/180;
 n = 100;
 [rho,theta] = meshgrid(vessel_profile_coords(1,:),linspace(vessel_sector_ini,vessel_sector_end,n));
 Z = repmat(vessel_profile_coords(2,:),n,1);
 X = rho.*cos(theta);
 Y = rho.*sin(theta);

 hold on
% plot3([0 0],[0 0],[-.2 .2],'b')
s_plot=surf(X,Y,Z);

% plot3(Trayectoria(1,200:end-1000),Trayectoria(2,200:end-1000),Trayectoria(3,200:end-1000),'g','linewidth',2)
plot3(Trayectoria(1,:),Trayectoria(2,:),Trayectoria(3,:),'g','linewidth',2)

% for i=1:N_toroidal_coils
% plot3(toroidal_coil_data(i).coords_coil(:,1),toroidal_coil_data(i).coords_coil(:,2)...
%     ,toroidal_coil_data(i).coords_coil(:,3),'r','linewidth',2)
% end
% plot3(toroidal_coil_data(1).coords_coil(:,1),toroidal_coil_data(1).coords_coil(:,2)...
%     ,toroidal_coil_data(1).coords_coil(:,3),'r')
% plot3(plasma_current_data(1).coords_coil(:,1),plasma_current_data(1).coords_coil(:,2)...
%     ,plasma_current_data(1).coords_coil(:,3),'c')
% plot3(lower_BvL_coil_data(1).coords_coil(:,1),lower_BvL_coil_data(1).coords_coil(:,2)...
%     ,lower_BvL_coil_data(1).coords_coil(:,3),'K')
% plot3(upper_BvL_coil_data(1).coords_coil(:,1),upper_BvL_coil_data(1).coords_coil(:,2)...
%     ,upper_BvL_coil_data(1).coords_coil(:,3),'K')
shading interp
s_plot.FaceColor=[1 1 1]*0.9;

camlight left

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
%title('ST25-D Particle Trajectory')
view(-14,35)
% view(0,90)
% view(0,0)
% view(29,0)
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;

axis([-1,1,-1,1,-.8,.8]*0.5)
box off
grid off
hold off
% s_plot.FaceAlpha=0.5;
%}


%%

% Plot poloidal cross section banana+Eq magnetico
%
% figure(2)
figure('pos',[600 0 500 500*1.1/.55])
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
    ,100*ones(length(toroidal_coil_data(1).coords_coil(:,2))),'k')
% plot3(plasma_center_xcoord,0,10,'xk','LineWidth',1.5)
% plot3(R_upper_BvL_coil,h_upper_BvL_coil,10,'xk','LineWidth',1.5)
% plot3(R_lower_BvL_coil,h_lower_BvL_coil,10,'xk','LineWidth',1.5)
% plot3(R_upper_BvU_coil,h_upper_BvU_coil,10,'xk','LineWidth',1.5)
% plot3(R_lower_BvU_coil,h_lower_BvU_coil,10,'xk','LineWidth',1.5)
plot3(vessel_profile_coords(1,:),vessel_profile_coords(2,:),10*ones(length(vessel_profile_coords)),'-k','LineWidth',1.5)
% plot3(plasma_profile(:,1),plasma_profile(:,2),10*ones(length(plasma_profile)),'--r','LineWidth',1.5)

phi_traject=atan2(Trayectoria(2,:),Trayectoria(1,:));

traject_2_plot_x=[Trayectoria(1,1:2:end)'./cos(phi_traject(1:2:end))'];
traject_2_plot_y=[Trayectoria(3,1:2:end)'];
traject_2_plot_z=100*ones(length(traject_2_plot_x),1);

plot3(traject_2_plot_x(700:end),traject_2_plot_y(700:end),traject_2_plot_z(700:end),'g','linewidth',1.5)



axis([0,.6,-.7,.7])
xlabel('x (m)')
ylabel('y (m)')

c = colorbar;
set(c,'fontsize',14);
ax = gca;
ax.FontSize = 14;
ylabel(c, 'Field modulus (T)');
caxis([0,0.7])
shading interp;
lighting gouraud;
% print('banana_20keV_2','-djpeg','-r900');
hold off
%}







