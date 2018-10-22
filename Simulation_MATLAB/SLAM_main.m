
% This is the main steps for SLAM
close all;
clear all;

num_steps = 10; %% number of time steps. delta t is assumed to be one.
dt = 1.;
% 1. we place the robot at one point, and that point is the origin of the global
% coordinate. I will assume the robot orientation is pi/2 with some noise.
% I will use robot_true as the true robot's x,y,theta,
% and use robot_initial as the guess initial.
robot_initial = [5/2*1000 0 pi/2]';
%robot_true = [3/2*1000 0 88/180*3.14]';
robot_true= robot_initial;
Cov_initial = diag([20^2 20^2 (2/180*3.14)^2]); %+-2degree std.

%%%% set the initial state
%x = robot_initial;
true_state = zeros(3, num_steps);
predicted_state = zeros(3, num_steps);
pure_propagation_state = zeros(3, num_steps);
error = zeros(2, num_steps);

%%%% Note: the unit used is mm, not m
%%%% set up the robot parameter
r = 250;
b = 300;
state_all = robot_initial;
Cov = Cov_initial;
endpoints_in_state = [];
x_true = robot_true;
for t = 1:num_steps
    
    if t ==1
        %sigma_u_x = 0;
        %sigma_u_y = 0;
        %sigma_u_theta = 0/180*3.14;
        sigma_v = 0;
        sigma_w = 0/180*3.14;
    else
        %sigma_u_x = 50;
        %sigma_u_y = 50;
        %sigma_u_theta = 5/180*3.14;
        sigma_v = 50;
        sigma_w = 5/180*3.14;
    end

    
    
    %%%% 1.Generate noise and record the true location
    %u_noise = [sigma_u_x sigma_u_y sigma_u_theta]'.*(2*rand(3,1)-1); % use uniform instead of gaussian to make noise bounded
    u_noise = [sigma_v sigma_w]'.*(2*rand(2,1)-1);
    u_input = [500 0]';
    u_true = u_input + u_noise;
    x_true = motion_model(x_true,u_true,dt);
    
    %% 2.apply the control input and predict your new state
    R_vw = diag([(2*sigma_v)^2/10,(2*sigma_w)^2/10]);
    [state_all,Cov_k_k_minus] = state_propagation(state_all,Cov,R_vw,u_input,dt);
    pure_propagation_state(:,t) = state_all(1:3);
    Cov = Cov_k_k_minus;
    
    %% 3.apply the laser scanner
    p_r = x_true(1:2,1)/1000;
    heading_r = x_true(3);
    sigma_d = 0.005;
    %sigma_d = 0;
    [d_m,theta_m] = gen_sim_laser_data(p_r,heading_r,sigma_d);
    theta_m = theta_m';
    d_m = d_m*1000;%change to mm unit
    Cov_single_point = [(sigma_d*1000)^2,0;0,1e-8];
    %% 4.Line Extraction
    % Note that the returned endpoints_in_line is in the line coordinate
    % not the global coordinate
    [line_para, endpoints_in_line,R_mat] = LineExtraction_main(d_m,theta_m,Cov_single_point); 
    
    %manually see which line corresponds to which
    %end_point_gl_m = [];
    for line_iter = 1:size(line_para,2)
        [end_point_1,end_point_2,end_point_gl] = end_point_2_global(state_all(1:3,1),endpoints_in_line(1:2,line_iter),endpoints_in_line(3:4,line_iter));
        plot([end_point_gl(1),end_point_gl(2)]'/1000,[end_point_gl(3),end_point_gl(4)]'/1000,'b-','LineWidth',5); 
    end
    
    %% 5.Line Matching
    if length(state_all)>3
       % only do line matching when state has lines
       state_line =state_all(4:end,1);
       line_in_state = reshape(state_line,[2,length(state_line)/2]);
       endpoints_in_state_use = endpoints_in_state;
       line_test = line_para;
       endpoints_test = endpoints_in_line;
       R_test = R_mat;
       robot_vector = state_all(1:3,1);
       orien_thresh = 22.5/180*3.14; %+-22.5deg
       dist_thresh = 1000;
       [line_robot,R_line,line_idx,endpoints_line_measure,line_un_match,endpoints_un_match,R_un_match] = LineMatch(line_in_state, endpoints_in_state_use,...
           line_test,endpoints_test,R_test,robot_vector,orien_thresh,dist_thresh);
       %% do the measurement update
       [state_all_new,Cov_k_k_new] = measurement_update(state_all,Cov,line_robot,R_line,line_idx);
       Cov = Cov_k_k_new;
       state_all = state_all_new;
       %% based on updated line parameters in state vector, update the matched lines ending points
       endpoints_in_state_new = update_mached_line_endpoint(state_all,endpoints_in_state,line_idx,endpoints_line_measure);
       endpoints_in_state = endpoints_in_state_new;
    end
    
    
    %% 6.add the lines that are not in the state vector
    if length(state_all)==3 || length(line_un_match)>0
        if length(state_all)==3
            line_add = line_para;
            R_add = R_mat;
            endpoints_add = endpoints_in_line;
        elseif length(line_un_match)>0
            line_add = line_un_match;
            R_add = R_un_match;
            endpoints_add = endpoints_un_match;
        end
            
        for i = 1:size(line_add,2)
            
             % Note we need to change the line parameters and endpoints to
            % the global coordinate before adding them
            % In order to change the line parameters to global coordinate,
            % we first change the ending points, which will then give us
            % the line parameter. Note, should not use single point
            % correspondence to determine the global coordinate.
            R_add_single = R_add(2*i-1:2*i,:);
            %line_add_theta = line_para(2,i);
            %line_add_d = line_para(1,i);
            %robot_vector = state_all(1:3,1);
            x_part = endpoints_add(1:2,i);
            y_part = endpoints_add(3:4,i);
            [state_all_new,Cov_k_k_new,end_point_gl] = add_line_2_state(state_all,Cov,x_part,y_part,R_add_single);

            plot([end_point_gl(1),end_point_gl(2)]'/1000,[end_point_gl(3),end_point_gl(4)]'/1000,'b-','LineWidth',8);  
            endpoints_in_state = [endpoints_in_state end_point_gl];
            state_all = state_all_new;
            Cov = Cov_k_k_new;
        end
    end
    predicted_state(:, t) = state_all(1:3,1);
    true_state(:, t) = x_true;
    
end
plot_endpoints_state(endpoints_in_state)        
        
       
