function [state_all_new,Cov_k_k_new,end_point_gl] = add_line_2_state(state_all,Cov,x_part,y_part,R_add_single)

robot_vector = state_all(1:3,1);
[~,~,end_point_gl] = end_point_2_global(robot_vector,x_part,y_part);

[end_point_1_r,end_point_2_r] = end_point_2_robot(x_part,y_part);
d_est = [norm(end_point_1_r),norm(end_point_2_r)]';
theta_est = [atan2(end_point_1_r(2),end_point_1_r(1)),atan2(end_point_2_r(2),end_point_2_r(1))]';
[d_l, theta_l] = line_para_cal(d_est,theta_est);
line_add_single_r = [d_l,theta_l]';
[state_all_new,Cov_k_k_new] = cov_update_new_line_measure(Cov,R_add_single,state_all,line_add_single_r);