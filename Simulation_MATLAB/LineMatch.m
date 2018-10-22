
function [line_robot,R_line,line_idx,endpoints_line_measure,line_un_match,endpoints_un_match,R_un_match] = LineMatch(line_in_state, endpoints_in_state,...
    line_test,endpoints_test,R_test,robot_vector,orien_thresh,dist_thresh)

% this function is to compare each pair's possible match based on 3
% criteria: orientation difference in polar coordinate, distance of
% mid-point of test line to the reference line, and the line overlapping
% It will split the newly obtained lines into 2 categories: matched,
% un-matched. 
% For matched lines, output the line parameters, R matrix, and the matched line position in state
% aliagned with the matched lines from the state. 
% Note: We should not update matched lines ending points!
% For un-matched lines, output the line parameters, R matrix, and ending
% points.



% Input: line_in_state, 2xN matrix, reference line or line already in state,
% each column represents one line configuration: d,theta
%        endpoints_in_state, 4xN matrix, reference line's ending points,
% each column represents 2 end points with the order x1,x2,y1,y2
%        line_test, 2xN matrix, in line coordinate. test line or newly obtained line with the
% laser. Each column represents one line configuration: d,theta.
%        endpoints_test, 4xN matrix, test line's ending points,
% each column represents 2 end points with the order x1,x2,y1,y2
%        R_test:2Nx2 matrix, each 2x2 matrix represents the covariance of
% the newly obtained estimated line parameters,d_star(the first one), theta_star
%        theta_turn, scalar in radian unit, estimated robot orientation
% w.r.t. the global coordinate
%        trans_move, 2x1 vector, robot translation w.r.t global coordinate
%        orien_thresh, scalar, radian unit, criteria to decide whether the
% two lines are the same line in terms of the orientation difference
%        dist_thresh, scalar, criteria to decide whether the two line are
%        same in terms of the distance

% Output: 
%        line_robot: 2xm matrix, each column is a matched measurement line parameter vector
%                    [d,theta]' in current robot coordinate. 
%       R_line: 2mx2 matrix, each 2x2 matrix is the covariance of the
%               line_robot parameters.
%       line_idx: mx1 vector, each element is the current line starting
%                 index
%        endpoints_line_measure, 4xN matrix, matched measurement line's
%                                ending points in line coordinate. each column represents 2 end points with the order x1,x2,y1,y2
%        line_un_match: 2xl matrix, un-matched lines . each column
%        represents one line configuration: d,theta
%        endpoints_un_match,4xl matrix, un-matched line's ending points,
%                            each column represents 2 end points with the order x1,x2,y1,y2
%        R_un_match: 2lx2 matrix, un-matched lines' covariance matrix.

trans_move = robot_vector(1:2,1);
theta_turn = robot_vector(3);
num_line_state = size(line_in_state,2);
num_line_test = size(line_test,2);

%% 1. change the test Line Parameters to robot coordinate and change the test line ending points to global coordinate 

end1_test_gl = [];% ending point 1
end2_test_gl = [];% ending point 2

line_test_r = [];

for k = 1:num_line_test
    x_part = endpoints_test(1:2,k);
    y_part = endpoints_test(3:4,k);
    %% 1.1 change the test line ending points 
    %robot_vector = trans_move;
    [end_point_1,end_point_2,~] = end_point_2_global(robot_vector,x_part,y_part);
    end1_test_gl = [end1_test_gl,end_point_1];
    end2_test_gl = [end2_test_gl,end_point_2];
    %% 1.2 change the test line parameters to robot coordinate
    [end_point_1_r,end_point_2_r] = end_point_2_robot(x_part,y_part);
    d_est = [norm(end_point_1_r),norm(end_point_2_r)]';
    theta_est = [atan2(end_point_1_r(2),end_point_1_r(1)),atan2(end_point_2_r(2),end_point_2_r(1))]';
    [d_l, theta_l] = line_para_cal(d_est,theta_est);
    line_add_single_r = [d_l,theta_l]';
    line_test_r = [line_test_r,line_add_single_r];
end



%% 2.want to compute each possible pair's fitting value
% reference is the row, test is the column
pair_val = zeros(num_line_state,num_line_test);
for i = 1:num_line_state
    % reference line already in the state
    line_single = line_in_state(:,i);
    theta_ref = line_single(2);
    d_ref = line_single(1);
    % note the turns are: x1, x2, y1, y2
    endpoint_ref = endpoints_in_state(:,i);
    
    for j = 1:num_line_test
        
        R_single = R_test(2*j-1:2*j,:);
        d_est_var = R_single(1,1);
        theta_est_var = R_single(2,2);
        theta_test = line_test_r(2,j);% in robot coordinate
        %theta_test = line_test_single(2);
        %d_test = line_test_single(1);

        
        %% 2.1 orientation difference
        
        % there are two possible measurement equations, select the lower
        % one. Note: need to scale the result to be between -pi and pi 
        dif_orien = [theta_test-theta_ref+theta_turn-pi;theta_test-theta_ref+theta_turn];
        
        % calculate the measurement equation for the distance to see which
        % measurement equation should be used
        d_test_m = [-d_ref+trans_move(1)*cos(theta_ref)+trans_move(2)*sin(theta_ref);d_ref-trans_move(1)*cos(theta_ref)-trans_move(2)*sin(theta_ref)];
        
        [~,measure_eq_idx] = max(d_test_m);
        dif_orien = abs(atan2(sin(dif_orien),cos(dif_orien)));
        
        
        %[dif_orien_use,measure_func_use] = min(abs(atan2(sin(dif_orien),cos(dif_orien))));
        %orien_com_val = dif_orien_use^2/theta_est_var;
        orien_com_val = dif_orien(measure_eq_idx);
        
        orien_val_sq_var = orien_com_val^2/theta_est_var;
        
        
        % normal distance difference
        %dif_d = abs([d_test+d_ref-trans_move(1)*cos(theta_ref)+trans_move(2)*sin(theta_ref);d_test-d_ref+trans_move(1)*cos(theta_ref)+trans_move(2)*sin(theta_ref)]);
        %dif_d_use = dif_d(measure_func_use);
        
        
        
        %% 2.2 compute the test line midpoint distance to the reference line
        end1_test_gl_single = end1_test_gl(:,j);
        end2_test_gl_single = end2_test_gl(:,j);
        mid_test_gl = end2_test_gl_single+0.5*(end1_test_gl_single-end2_test_gl_single);
        %mid_test_gl = Rot_rob_12*(mid_test_rob-trans_move);
        dist_com_val = abs(cos(theta_ref)*mid_test_gl(1)+sin(theta_ref)*mid_test_gl(2)-d_ref);
        
        %dist_val_sq_var = dist_com_val^2/d_est_var;
  
        %% 2.3 compute the overlapping rate
        
        % use the current estimated orientation and translation of robot to
        % change the test points(which is w.r.t. current robot pos) to the
        % coordinate w.r.t. global coordinate. And then project them back
        % to the reference line to decide the overlapping rate
        
        %end1_test_gl = Rot_rob_12*(end1_test_rob-trans_move);
        %end2_test_gl = Rot_rob_12*(end2_test_rob-trans_move);
        
        end1_test_gl_in_line = zeros(2,1);
        end1_test_gl_in_line(1,1) = sin(theta_ref)*(sin(theta_ref)*end1_test_gl_single(1,1)-cos(theta_ref)*end1_test_gl_single(2,1))...
            -cos(theta_ref)*(-1*d_ref);
        end1_test_gl_in_line(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*end1_test_gl_single(1,1)+cos(theta_ref)*end1_test_gl_single(2,1))...
            -sin(theta_ref)*(-1*d_ref);
        
        % use the relative angle w.r.t. the line parameter theta to
        % determine the relative position of the 4 points 
        end1_test_theta = atan2(end1_test_gl_in_line(2,1),end1_test_gl_in_line(1,1))-theta_ref;
        end1_test_theta = atan2(sin(end1_test_theta),cos(end1_test_theta));
            
        end2_test_gl_in_line = zeros(2,1);
        end2_test_gl_in_line(1,1) = sin(theta_ref)*(sin(theta_ref)*end2_test_gl_single(1,1)-cos(theta_ref)*end2_test_gl_single(2,1))...
            -cos(theta_ref)*(-1*d_ref);
        end2_test_gl_in_line(2,1) = cos(theta_ref)*(-1*sin(theta_ref)*end2_test_gl_single(1,1)+cos(theta_ref)*end2_test_gl_single(2,1))...
            -sin(theta_ref)*(-1*d_ref);
        
        end2_test_theta = atan2(end2_test_gl_in_line(2,1),end2_test_gl_in_line(1,1))-theta_ref;  
        end2_test_theta = atan2(sin(end2_test_theta),cos(end2_test_theta));
        %if end2_test_theta <0
        %    end2_test_theta = end2_test_theta + pi;
        %end
        
        end1_ref_theta = atan2(endpoint_ref(3),endpoint_ref(1))-theta_ref;
        end1_ref_theta = atan2(sin(end1_ref_theta),cos(end1_ref_theta));
        %if end1_ref_theta<0
        %    end1_ref_theta = end1_ref_theta + pi;
        %end
        
        end2_ref_theta = atan2(endpoint_ref(4),endpoint_ref(2))-theta_ref;
        end2_ref_theta = atan2(sin(end2_ref_theta),cos(end2_ref_theta));
        %if end2_ref_theta<0
        %    end2_ref_theta = end2_ref_theta+pi;
        %end
        
        
        if end1_test_theta > end2_test_theta
            temp = end1_test_theta;
            end1_test_theta = end2_test_theta;
            end2_test_theta = temp;
        end
        
        if end1_ref_theta > end2_ref_theta
            temp = end1_ref_theta;
            end1_ref_theta = end2_ref_theta;
            end2_ref_theta = temp;
        end
        
        if end1_test_theta >= end1_ref_theta && end2_test_theta<=end2_ref_theta % test line segment within the reference segment
            lap_rate = 1;
        elseif end1_test_theta>=end1_ref_theta && end1_test_theta<=end2_ref_theta % test seg starting point within reference seg, while end not 
            lap_rate = (end2_ref_theta-end1_test_theta)/(end2_ref_theta-end1_ref_theta);
        elseif end2_test_theta>=end1_ref_theta && end2_test_theta<=end2_ref_theta % test seg end point within reference seg, while start not
            lap_rate = (end2_test_theta-end1_ref_theta)/(end2_ref_theta-end1_ref_theta);
        elseif end1_ref_theta>= end1_test_theta && end2_ref_theta<= end2_test_theta && end2_ref_theta >= end1_test_theta % test seg encompass the reference seg
            lap_rate = 1;
        else
            lap_rate = 0;
        end
        
        if orien_com_val<=orien_thresh &&dist_com_val<=dist_thresh&& lap_rate>=0.05
            pair_val(i,j) = 0.2*orien_com_val/orien_thresh+0.3*dist_com_val/dist_thresh+0.5*(1-lap_rate)/0.95;
        else
            pair_val(i,j) = 10;% set a large value, valid value should be <=1
        end
        % use the variance-weighted difference as criteria
%         if orien_val_sq_var<=orien_thresh &&dist_val_sq_var<=dist_thresh&& lap_rate>=0.2
%             pair_val(i,j) = 0.2*orien_val_sq_var/orien_thresh+0.3*dist_val_sq_var/dist_thresh+0.5*(1-lap_rate)/0.8;
%         else
%             pair_val(i,j) = 10;% set a large value, valid value should be <=1
%         end


    end
    
end

%% 3. obtain the matching based on the pair_val
pair_match = [];
for i = 1:num_line_state
    pair_val_single = pair_val(i,:);
    if min(pair_val_single)<=1
        [~,idx_match] = min(pair_val_single);
        pair_match = [pair_match;i, idx_match];
    end
end

%% 4. update the ending points for the matched lines in state
line_idx = zeros(size(pair_match,1),1); % line idx corresponding to the function measurement_update
%line_gl = zeros(2,size(pair_match,1)); % 2xm line parameters of the matched ones from state
line_robot = zeros(2,size(pair_match,1)); %2xm matched lines from the test(new measurement)
R_line = zeros(2*size(pair_match,1),2); % 2mx2 covariance matrix for the matched lines from the test
endpoints_line_measure = zeros(4,size(pair_match,1));
for i = 1:size(pair_match,1)
    
    ref_idx = pair_match(i,1);
    test_idx = pair_match(i,2);
 
    line_idx(i,1) = (ref_idx-1)*2+1+3; % add 3 in order to take into the robot 3 parameters into account
    %line_gl(1,i) = d_ref;
    %line_gl(2,i) = theta_ref;
    line_robot(:,i) = line_test_r(:,test_idx);
    R_line(2*i-1:2*i,:) = R_test(2*test_idx-1:2*test_idx,:);
    endpoints_line_measure(:,i) = endpoints_test(:,test_idx);
end

%% 5 get the un-matched lines
line_un_match = line_test;
idx_match = unique(pair_match(:,2)); % unique idx sorted from smallest to largest
line_un_match(:,idx_match) = [];
endpoints_un_match = endpoints_test;
endpoints_un_match(:,idx_match) = [];
R_un_match = R_test;
for i = length(idx_match):-1:1
    R_un_match(2*idx_match(i)-1:2*idx_match,:) = [];
end





        