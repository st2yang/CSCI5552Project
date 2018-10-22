function [line_para, endpoints_in_line,R_mat] = LineExtraction_main(d_est,theta_est,Cov_Laser)
% This function is the main call function to do the Line Extraction
% its outputs are the line parameter, 
% the end points on the line(not the end points that are usually not on the line)

% Input: d_est, nx1 vector, the distance measurements for each point
%        theta_est, nx1 vector, the angle measurements for each point
%        Cov_Laser:2x2,covariance of the single laser point,[d,theta] order

% Output: line_para,2xN matrix, each column represents 2 parameters (d,theta) of a line in
%                   2D space (polar coordinate).
%         endpoints_in_line,4xN matrix, each column represents end points
%                           in the order of (x1,x2,y1,y2)

%         R_mat, 2Nx2 matrix, each 2x2 matrix represents the covariance of
%                the estimated line parameters,d_star(the first one), theta_star
% 

maxIter = 1000;
min_num = 10;
thresh = 250; %maximum distance from point to the line
min_rest = 20;
theta_thresh = 1;
seg_dist_thresh = 500;
line_seg_dist_thresh = 500;
% 1. RANSAC
%[line_RANSAC,endpoints,R_mat] = LineExtraction_RANSAC(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh);

[line_RANSAC,endpoints,R_mat] = LineExtraction_RANSAC_2(d_est,theta_est,thresh,min_num,maxIter,min_rest,theta_thresh,seg_dist_thresh,Cov_Laser);


%[line_RANSAC,endpoints,R_mat,~] = LineSegment_Diff_method(d_est,theta_est,seg_dist_thresh,min_num);

endpoints_in_line = zeros(size(endpoints)); % first 2 are x_axis, next 2 are y_axis
idx_remove = [];
for i = 1:size(line_RANSAC,2)
    %obtain each line parameters
    cos_theta = cos(line_RANSAC(2,i));
    sin_theta = sin(line_RANSAC(2,i));
    d_theta = line_RANSAC(1,i);
    %project the start and end points to the line
    start_x = endpoints(1,i)*cos(endpoints(2,i));
    start_y = endpoints(1,i)*sin(endpoints(2,i));
    end_x = endpoints(3,i)*cos(endpoints(4,i));
    end_y = endpoints(3,i)*sin(endpoints(4,i));
    endpoints_in_line(1,i) = sin_theta*(sin_theta*start_x-cos_theta*start_y)-cos_theta*(-1*d_theta);
    endpoints_in_line(3,i) = cos_theta*(-1*sin_theta*start_x+cos_theta*start_y)-sin_theta*1*(-1*d_theta);
    
    endpoints_in_line(2,i) = sin_theta*(sin_theta*end_x-cos_theta*end_y)-cos_theta*(-1*d_theta);
    endpoints_in_line(4,i) = cos_theta*(-1*sin_theta*end_x+cos_theta*end_y)-sin_theta*(-1*d_theta);
    if sqrt((endpoints_in_line(1,i)-endpoints_in_line(2,i))^2+(endpoints_in_line(3,i)-endpoints_in_line(4,i))^2)<line_seg_dist_thresh
        idx_remove = [idx_remove;i];
    end
end
if length(idx_remove)>0
    line_RANSAC(:,idx_remove) = [];
    endpoints_in_line(:,idx_remove) = [];
    for i = length(idx_remove):-1:1 % remove the highest index first to make sure the index is correct
        R_mat(2*idx_remove(i)-1:2*idx_remove(i),:) = [];
    end
end

line_para = line_RANSAC;
