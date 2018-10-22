function  [d_H, theta_H] = line_para_cal(d_est,theta_est)
% compute the line parameters by least square 

u_stack = [d_est.*cos(theta_est), d_est.*sin(theta_est)];
u_vec = mean(u_stack,1)'; %2x1 vector
dif_u_p = u_vec*ones(1,length(theta_est))-u_stack';
A = dif_u_p*dif_u_p';
[V,D] = eig(A);
[~,idx_min] = min(diag(D));
V_min = V(:,idx_min);
cos_theta= V_min(1);
sin_theta= V_min(2);
d_theta = u_vec'*V_min;
coef_est = [cos_theta,sin_theta,d_theta]';
coef_est = sign(d_theta)*coef_est;
theta_H = atan2(coef_est(2),coef_est(1));
d_H = coef_est(3);