%% TEST: D-shaped Coil

z2=0.896/2;
R1=0.7;
n_div_arc1=10;
R2=0.12;
n_div_arc2=5;
x_tang=.2;
center_R1_x_coord=-.25;
x_max=R1+center_R1_x_coord;

arc1_x_coord_up=x_tang+.02:(x_max-x_tang-.02)/n_div_arc1:x_max;
arc1_y_coord_up=sqrt(R1^2-(arc1_x_coord_up-center_R1_x_coord).^2);
% arc1_x_coord=[arc1_x_coord_up arc1_x_coord_up(end:-1:1)];
% arc1_y_coord=[arc1_y_coord_up -arc1_y_coord_up(end:-1:1)];


arc2_x_coord_up=0:x_tang/n_div_arc2:x_tang;
arc2_y_coord_up=z2+sqrt(R2^2-(arc2_x_coord_up-R2).^2);
% arc2_x_coord=[arc2_x_coord_up arc2_x_coord_up(end:-1:1)];
% arc2_y_coord=[arc2_y_coord_up -arc2_y_coord_up(end:-1:1)];


Dcoil_profile_x_coords=[arc2_x_coord_up arc1_x_coord_up arc1_x_coord_up(end:-1:1) arc2_x_coord_up(end:-1:1) arc2_x_coord_up(1)];
Dcoil_profile_y_coords=[arc2_y_coord_up arc1_y_coord_up -arc1_y_coord_up(end:-1:1) -arc2_y_coord_up(end:-1:1) arc2_y_coord_up(1)];

% plot(Dcoil_profile_x_coords,Dcoil_profile_y_coords,'-r')
% axis equal