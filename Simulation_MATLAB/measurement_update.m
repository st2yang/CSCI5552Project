
function [state_all_new,Cov_k_k_new] = measurement_update(state_all,Cov_k_k_minus,line_robot,R_line,line_idx)

% In this function, measurement update is performed according the line matching
% results.

% Input:
%       state_all: (3+2n)x1 state vector, first 3 are robot (x,y,theta)'

%       Cov_k_k_minus: covariance matrix for the state_all vector

%       line_robot: 2xm matrix, each column is a line parameter vector
%       [d,theta]' in current robot coordinate. And it is already aliagned
%       with the line_idx vector to current state vector.

%       R_line: 2mx2 matrix, each 2x2 matrix is the covariance of the
%       line_robot parameters.

%       line_idx: mx1 vector, each element is the current line starting
%       index

% Output:
%        state_all_new: (3+2n)x1 state vector after measurement update

%        Cov_k_k_new: (3+2n)x(3+2n) covariance matrix after measurement
%        update
H_all = [];
measure_dif_all = [];
num_state = length(state_all);
R_all = [];
for i = 1:length(line_idx)
    robot_gl = state_all(1:3);
    %line_gl_single = line_gl(:,i);
    line_gl_single = state_all(line_idx(i):line_idx(i)+1,1);
    line_idx_single = line_idx(i);
    [H_single,eq_idx] = cal_H_mat(robot_gl,line_gl_single,line_idx_single,num_state);
    line_robot_single = line_robot(:,i);
    if eq_idx ==1
        orien_dif = line_robot_single(2)-line_gl_single(2)+robot_gl(3);
        measure_dif = [line_robot_single(1)-line_gl_single(1)+robot_gl(1)*cos(line_gl_single(2))+robot_gl(2)*sin(line_gl_single(2));...
            atan2(sin(orien_dif),cos(orien_dif))];
    else
        orien_dif = line_robot_single(2)-line_gl_single(2)+robot_gl(3)-pi;
        measure_dif = [line_robot_single(1)+line_gl_single(1)-robot_gl(1)*cos(line_gl_single(2))-robot_gl(2)*sin(line_gl_single(2));...
            atan2(sin(orien_dif),cos(orien_dif))];
    end
    measure_dif_all = [measure_dif_all;measure_dif];
    H_all = [H_all;H_single];
    R_all = blkdiag(R_all,R_line(2*i-1:2*i,:));
end
% do the EKF measurement update
S = H_all*Cov_k_k_minus*H_all'+R_all;
K = Cov_k_k_minus*H_all'*inv(S);
% update the covariance
Cov_k_k_new = Cov_k_k_minus-K*H_all*Cov_k_k_minus;
% enforce symmetric
Cov_k_k_new = (Cov_k_k_new+Cov_k_k_new')/2;
% update the estimated variables
state_all_new = state_all+K*measure_dif_all;
    
    
    
    