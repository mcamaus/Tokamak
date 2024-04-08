%Poincare plot at a given phi: Execute field_line_tracing before

%
%
phi_poincare=0*pi/180;

for i=1:n_lines_z
    for j=1:n_lines_x
        traject(i,j).phi=atan2(traject(i,j).coords(:,2),traject(i,j).coords(:,1));
        %set the index of poincare points to 1
        np=1;
        for k=1:length(traject(i,j).phi)-1
            % detect the cros-section of the trajectory with the plane y1-y2
            if(traject(i,j).phi(k)<=phi_poincare && traject(i,j).phi(k+1)>phi_poincare)
                % store detected cross-section point
                traject(i,j).ps_coords2(np,:)=traject(i,j).coords(k,:);
                [traject(i,j).ps_coords(np,:),check]=plane_line_intersect(...
                    [-sin(phi_poincare) cos(phi_poincare) 0],[0 0 0],...
                    traject(i,j).coords(k,:),traject(i,j).coords(k+1,:));
                
                %         increase the index of poincare point
                np=np+1;
            end
        end
    end
end

% Filter to plot:
% Mesh filter in position


% Filter n� of loops

%{
% clf
figure('pos',[300 100 450 800])
hold on
plot3(toroidal_coil_data(1).coords_coil(:,1),toroidal_coil_data(1).coords_coil(:,3)...
    ,100.*10*ones(length(toroidal_coil_data(1).coords_coil(:,2))),'k')
% plot3(plasma_center_xcoord,0,10,'xk','LineWidth',1.5)
% plot3([x_min x_min x_max x_max x_min],[z_min z_max z_max z_min z_min],[0 0 0 0 0],'--g','LineWidth',1)
plot3(R_upper_BvL_coil,h_upper_BvL_coil,10,'xk','LineWidth',1.5)
plot3(R_lower_BvL_coil,h_lower_BvL_coil,10,'xk','LineWidth',1.5)
plot3(R_upper_BvU_coil,h_upper_BvU_coil,10,'xk','LineWidth',1.5)
plot3(R_lower_BvU_coil,h_lower_BvU_coil,10,'xk','LineWidth',1.5)
plot3(R_upp_DIV_coil,h_upp_DIV_coil,10,'xk','LineWidth',1.5)
plot3(R_low_DIV_coil,h_low_DIV_coil,10,'xk','LineWidth',1.5)

plot3(vessel_profile_coords(1,:),vessel_profile_coords(2,:),10*ones(length(vessel_profile_coords)),'-k','LineWidth',1.5)
% plot3(plasma_profile(:,1),plasma_profile(:,2),10*ones(length(plasma_profile)),'--r','LineWidth',1)
for i=1:n_lines_z
    for j=1:n_lines_x
        if isempty(traject(i,j).ps_coords) ==1
            
%         elseif traject(i,j).ps_coords(:,3) > .3
            
        else
            plot( traject(i,j).ps_coords(:,1)/cos(phi_poincare), traject(i,j).ps_coords(:,3),'b.','markersize', 5);  %poincare points
        end
        
%         plot( traject(i,j).ps_coords(:,1), traject(i,j).ps_coords(:,3),'b.','markersize', 5);
%         plot( traject(i,j).coords(1,1), traject(i,j).coords(1,3),'rx','markersize', 5);  % initial poins of trajectories
        %     plot( traject(k).ps_coords2(:,1),traject(k).ps_coords2(:,3),'r.','markersize', 5);  % poincare points without intersection
    end
end
view(0,90);
axis([0,.55,-.5,.5])
xlabel('radial coordinate (m)');
ylabel('vertical coordinate (m)');
hold off
%}

%Field line plot
%
%Plot field lines
figure(2)
hold on
i_plot=2;
j_plot=2;
plot3([0 0],[0 0],[-.2 .2],'b')
plot3([0 .5*cos(phi_poincare)],[0 .5*sin(phi_poincare)],[0 0],'k') %poincare plane
if isempty(traject(i_plot,j_plot).ps_coords) ==1 
else
plot3(traject(i_plot,j_plot).ps_coords(:,1),...
    traject(i_plot,j_plot).ps_coords(:,2),...
    traject(i_plot,j_plot).ps_coords(:,3),'xr','markersize', 10) %poincare points
end
% plot3([traject(i_plot).ps_coords2(:,1)],...
%     [traject(i_plot).ps_coords2(:,2)],...
%     [traject(i_plot).ps_coords2(:,3)],'xc','markersize', 10) %poincare points2


plot3(traject(i_plot,j_plot).coords(:,1),traject(i_plot,j_plot).coords(:,2),traject(i_plot,j_plot).coords(:,3),'b','linewidth',2)

for i=1:N_toroidal_coils
plot3(toroidal_coil_data(i).coords_coil(:,1),toroidal_coil_data(i).coords_coil(:,2)...
    ,toroidal_coil_data(i).coords_coil(:,3),'r','linewidth',2)
end
% plot3(plasma_current_data(1).coords_coil(:,1),plasma_current_data(1).coords_coil(:,2)...
%     ,plasma_current_data(1).coords_coil(:,3),'c')
% plot3(lower_BvL_coil_data(1).coords_coil(:,1),lower_BvL_coil_data(1).coords_coil(:,2)...
%     ,lower_BvL_coil_data(1).coords_coil(:,3),'K')
% plot3(upper_BvL_coil_data(1).coords_coil(:,1),upper_BvL_coil_data(1).coords_coil(:,2)...
%     ,upper_BvL_coil_data(1).coords_coil(:,3),'K')

plot3(lower_BvL_coil_data(1).coords_coil(:,1),lower_BvL_coil_data(1).coords_coil(:,2)...
    ,lower_BvL_coil_data(1).coords_coil(:,3),'K','linewidth',2)
plot3(upper_BvL_coil_data(1).coords_coil(:,1),upper_BvL_coil_data(1).coords_coil(:,2)...
    ,upper_BvL_coil_data(1).coords_coil(:,3),'K','linewidth',2)
plot3(lower_BvU_coil_data(1).coords_coil(:,1),lower_BvU_coil_data(1).coords_coil(:,2)...
    ,lower_BvU_coil_data(1).coords_coil(:,3),'K','linewidth',2)
plot3(upper_BvU_coil_data(1).coords_coil(:,1),upper_BvU_coil_data(1).coords_coil(:,2)...
    ,upper_BvU_coil_data(1).coords_coil(:,3),'K','linewidth',2)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('ST25-D field line')
view(-14,4)
% view(0,90)
% view(0,0)

axis([-1,1,-1,1,-1,1]*0.7)
box
grid
hold off
%}


%Plot poincare + magnetic equilibrium
%
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
% plot3(plasma_center_xcoord,0,10,'xk','LineWidth',1.5)
% plot3(R_upper_BvL_coil,h_upper_BvL_coil,10,'xk','LineWidth',1.5)
% plot3(R_lower_BvL_coil,h_lower_BvL_coil,10,'xk','LineWidth',1.5)
% plot3(R_upper_BvU_coil,h_upper_BvU_coil,10,'xk','LineWidth',1.5)
% plot3(R_lower_BvU_coil,h_lower_BvU_coil,10,'xk','LineWidth',1.5)
% plot3(R_upp_DIV_coil,h_upp_DIV_coil,10,'xk','LineWidth',1.5)
% plot3(R_low_DIV_coil,h_low_DIV_coil,10,'xk','LineWidth',1.5)

plot3(vessel_profile_coords(1,:),vessel_profile_coords(2,:),10*ones(length(vessel_profile_coords)),'-k','LineWidth',1.5)
% plot3(plasma_profile(:,1),plasma_profile(:,2),10*ones(length(plasma_profile)),'--r','LineWidth',1.5)

for i=1:n_lines_z
    for j=1:n_lines_x
        if isempty(traject(i,j).ps_coords) ==1
            
        else
            plot3( traject(i,j).ps_coords(:,1)/cos(phi_poincare), traject(i,j).ps_coords(:,3),100*ones(length(traject(i,j).ps_coords(:,1))),'g.','markersize', 8);  %poincare points
        end
        
%         plot( traject(i,j).ps_coords(:,1), traject(i,j).ps_coords(:,3),'b.','markersize', 5);
%         plot( traject(i,j).coords(1,1), traject(i,j).coords(1,3),'rx','markersize', 5);  % initial poins of trajectories
        %     plot( traject(k).ps_coords2(:,1),traject(k).ps_coords2(:,3),'r.','markersize', 5);  % poincare points without intersection
    end
end

axis([0,.6,-.7,.7])
xlabel('radial coordinate (m)');
ylabel('vertical coordinate (m)');
c = colorbar;
set(c,'fontsize',14);
ax = gca;
ax.FontSize = 14;
ylabel(c, 'Field modulus (T)');
caxis([0,.8])
shading interp;
lighting gouraud;
% print('poincare_1','-dpng','-r900');
hold off
%}










