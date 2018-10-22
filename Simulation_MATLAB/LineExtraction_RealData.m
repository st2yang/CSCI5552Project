
clear all;
close all
%fileID = fopen('LaserData_forward.txt','r');
%fileID = fopen('LaserData_left.txt','r');
fileID = fopen('LaserData_right.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
dataNum = 361;
A = reshape(A,2,dataNum);
A = A';
x = A(:,1); %unit is mm
y = A(:,2); %unit is mm

plot(x,y,'r*','MarkerSize',10);hold on

% 2. translate points into polar coordinate
d_est = (x.^2 + y.^2).^0.5;
theta_est = atan2(y,x);

%% line extraction
Cov_single_point = [25,0;0,1e-4];%based on the material of 5551
[line_para_1, endpoints_in_line_1,R_mat_1] = LineExtraction_main(d_est,theta_est,Cov_single_point);
%plot the approximated line segments
for i = 1:size(line_para_1,2)
    x_part = endpoints_in_line_1(1:2,i)';
    y_part = endpoints_in_line_1(3:end,i)';
    plot(x_part,y_part,'b-');
    hold on;
end
%hold off;

%% make some transformation of the data points to mimic the move of the robot
%figure;
theta_move = 30/180*3.14+5/180*3.14*randn;%assume it turn 30deg with noise std 5 deg
trans_move = [0;1800] + 80*randn(2,1);%assume it move forward 1800mm with noise

%% make the points to be robot-coordinate
g_R_r = [cos(theta_move), -1*sin(theta_move) trans_move(1);sin(theta_move), cos(theta_move) trans_move(2); 0 0 1];
Point_pos_g = [x';y';ones(1,length(x))];
Point_pos_r = inv(g_R_r)*Point_pos_g;
x2 = Point_pos_r(1,:)';
y2 = Point_pos_r(2,:)';

% Rot_mat_12 = [cos(theta_move), -1*sin(theta_move);sin(theta_move), cos(theta_move)]; %This is 1R2, we want 2R1
% Rot_mat_21 = Rot_mat_12';
% Point_pos_1 = [x';y'];
% Point_pos_2 = Rot_mat_21*Point_pos_1+trans_move*ones(1,size(Point_pos_1,2));
% x2 = Point_pos_2(1,:)';
% y2 = Point_pos_2(2,:)';

% change the current points to the global coordinate(the previous coordinate)

% plot(x2,y2,'go','MarkerSize',10);
% hold on;

% translate points into polar coordinate
d_est2 = (x2.^2 + y2.^2).^0.5;
theta_est2 = atan2(y2,x2);

%% line extraction
% Note: the returned lines are in line-coordinate(in this simulation, it is actually robot-coordinate)
[line_para_2, endpoints_in_line_2, R_mat_2] = LineExtraction_main(d_est2,theta_est2,Cov_single_point);

%% change the lines' parameters and endpoints to the global coordinate

Rot_line_2_rob = [cos(-pi/2), -1*sin(-pi/2) 0;sin(-pi/2), cos(-pi/2), 0;0 0 1];
Rot_line_2_rob = eye(3);
g_R_r = [cos(theta_move), -1*sin(theta_move) trans_move(1);sin(theta_move), cos(theta_move) trans_move(2); 0 0 1];

line_para_2_gl = [];
endpoints_in_line_2_gl = [];
for i = 1:size(line_para_2,2)
    d_l = line_para_2(1,i);
    theta_l = line_para_2(2,i);
    l_g = g_R_r*Rot_line_2_rob*[d_l*cos(theta_l), d_l*sin(theta_l),1]';
    line_para_2_gl = [line_para_2_gl,[norm(l_g(1:2)),atan2(l_g(2),l_g(1))]'];
    end_p = zeros(4,1);
    for j = 1:2
        end1_line = [endpoints_in_line_2(j,i),endpoints_in_line_2(j+2,i),1]';
        end1_gl = g_R_r*Rot_line_2_rob*end1_line;
        end_p(j) = end1_gl(1);
        end_p(j+2) = end1_gl(2);
    end
    endpoints_in_line_2_gl = [endpoints_in_line_2_gl,end_p];
end
    

%plot the approximated line segments
for i = 1:size(line_para_2,2)
    x_part = endpoints_in_line_2_gl(1:2,i)';
    y_part = endpoints_in_line_2_gl(3:end,i)';
    plot(x_part,y_part,'b-','LineWidth',5);
    hold on;
end
hold off;

%% Line matching
line_para_base = line_para_1;
endpoints_base = endpoints_in_line_1;
theta_guess = 30/180*3.14;
trans_move_guess = [0;1800];

% threshold setting, initial guess, should be tuned in real experiment
orien_thresh = 22.5/180*3.14; %+-22.5deg
dist_thresh = 1000;
% orien_thresh = 5e3;
% dist_thresh = 1000;
pair = LineMatch(line_para_base,endpoints_base,line_para_2,endpoints_in_line_2,R_mat_2,theta_guess,trans_move_guess,orien_thresh,dist_thresh);


% maxIter = 1000;
% min_num = 6;
% thresh = 30; %maximum distance from point to the line
% min_rest = 6;
% theta_thresh = 3;
% seg_dist_thresh = 400;
% % 1. RANSAC
% [line_RANSAC,endpoints,R_mat] = LineExtraction_RANSAC(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh);
% endpoints_in_line = zeros(size(endpoints)); % first 2 are x_axis, next 2 are y_axis
% for i = 1:size(line_RANSAC,1)
%     %obtain each line parameters
%     cos_theta = line_RANSAC(i,1);
%     sin_theta = line_RANSAC(i,2);
%     d_theta = line_RANSAC(i,3);
%     %project the start and end points to the line
%     start_x = endpoints(1,i)*cos(endpoints(2,i));
%     start_y = endpoints(1,i)*sin(endpoints(2,i));
%     end_x = endpoints(3,i)*cos(endpoints(4,i));
%     end_y = endpoints(3,i)*sin(endpoints(4,i));
%     endpoints_in_line(1,i) = sin_theta*(sin_theta*start_x-cos_theta*start_y)-cos_theta*(-1*d_theta);
%     endpoints_in_line(3,i) = cos_theta*(-1*sin_theta*start_x+cos_theta*start_y)-sin_theta*1*(-1*d_theta);
%     
%     endpoints_in_line(2,i) = sin_theta*(sin_theta*end_x-cos_theta*end_y)-cos_theta*(-1*d_theta);
%     endpoints_in_line(4,i) = cos_theta*(-1*sin_theta*end_x+cos_theta*end_y)-sin_theta*(-1*d_theta);
%     
% end


    
    
    
    
    
