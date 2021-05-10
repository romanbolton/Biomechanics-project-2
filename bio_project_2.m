clear all
close all
load 'bending.mat' 
%load 'bending_test.mat'

%% restore old matlab colormap and hold function
addpath(genpath('fix hold'));% hold function is adapted (see hold function current folder) such that the color order index is not reset (new plot starts with blue again)
set(groot,'defaultAxesColorOrder',[0 0 1;0 0.5 0;1 0 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25])

figure;
plot (acc);
ylabel('Acceleration (m/s^2)'); xlabel('Timesample');
title('Bending: Acceleration of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);     %Change the name of the title for every motion
legend('A_T_r_u_n_k_x', 'A_T_r_u_n_k_y', 'A_T_r_u_n_k_z','fontsize',20);

corrected_acc = acc((201:1000),:);

figure;
plot(corrected_acc);
ylabel('Acceleration (m/s^2)'); xlabel('Timesample');
title('Bending: Acceleration of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);     %Change the name of the title for every motion
legend('A_T_r_u_n_k_x', 'A_T_r_u_n_k_y', 'A_T_r_u_n_k_z','fontsize',20)

G_ref = corrected_acc(1,:);
G_cur = corrected_acc;
G_dif = G_cur - G_ref;
 
L_G_ref = pytha_col(G_ref);
L_G_cur = pytha_col(G_cur);
L_G_dif = pytha_col(G_dif);

gamma_rad = acos ((L_G_ref.^2 + L_G_cur.^2 - L_G_dif.^2)./(2.*L_G_ref.*L_G_cur));

gamma_deg = rad2deg (gamma_rad);

figure
plot (gamma_deg)
ylabel('Angles (degrees)'); xlabel('Timesample');
title('Bending: Angle of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);     %Change the name of the title for every motion

%step 3.3

t_old = t(201:1000);
fs_intended = round(1/nanmedian(diff(t_old)));
t_new = t_old(1):1/fs_intended:t_old(end);
gap_limit = 1;
corrected_acc = spline_interp_find_gaps (corrected_acc, gap_limit, t_old, t_new);

figure
plot (corrected_acc)
ylabel('Acceleration (m/s^2)'); xlabel('Timesample');
title('Bending: Acceleration of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);     %Change the name of the title for every motion
legend('A_T_r_u_n_k_x', 'A_T_r_u_n_k_y', 'A_T_r_u_n_k_z','fontsize',20);

fs=50;
Fc=1;

[smooth] = filtcol_new(corrected_acc,fs,Fc,[],[],[],1);

figure
plot (smooth)
ylabel('Acceleration (m/s^2)'); xlabel('Timesample');
title('Bending: Acceleration of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);     %Change the name of the title for every motion
legend('A_T_r_u_n_k_x', 'A_T_r_u_n_k_y', 'A_T_r_u_n_k_z','fontsize',20);

G_ref_smooth = smooth(1,:);
G_cur_smooth = smooth;
G_dif_smooth = G_cur_smooth - G_ref_smooth;
 
L_G_ref_smooth = pytha_col(G_ref_smooth);
L_G_cur_smooth = pytha_col(G_cur_smooth);
L_G_dif_smooth = pytha_col(G_dif_smooth);

gamma_rad_smooth = acos ((L_G_ref_smooth.^2 + L_G_cur_smooth.^2 - L_G_dif_smooth.^2)./(2.*L_G_ref_smooth.*L_G_cur_smooth));

gamma_deg_smooth = rad2deg (gamma_rad_smooth);

figure
plot (gamma_deg_smooth)
ylabel('Angles (degrees)'); xlabel('Timesample');
title('Bending: Angle of the trunk'); set(gca, 'Box', 'off', 'fontsize',20);











