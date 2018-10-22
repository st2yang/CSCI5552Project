function [H,eq_idx] = cal_H_mat(robot_gl,line_gl,line_idx,num_state)

% This function calculate the H matrix used in EKF.
% Note this H is only one line's H matrix

% Input:
%        robot_gl: 3x1 vector, robot (x,y,theta)' in global coordinate
%        line_gl: 2x1 vector,line global coordinate parameters (d,theta)
%        line_idx: scalar, current line starting idx
%        num_state: scalar, current number of states

% Output:
%        H: 2x(3+2n) matrix, Jacobian matrix of the single measurement
%        eq_idx: scalar, measurement equation used

% 1.decide which measurement equation should be used
d_line = zeros(2,1);
d_line(1) = line_gl(1)-robot_gl(1)*cos(line_gl(2))-robot_gl(2)*sin(line_gl(2));
d_line(2) = -line_gl(1)+robot_gl(1)*cos(line_gl(2))+robot_gl(2)*sin(line_gl(2));
[~,eq_idx] = max(d_line);

% 2.calculate the H matrix w.r.t. 2 possible measurement equation
H = zeros(2,num_state);
if eq_idx==1
    H(1,1) = -cos(line_gl(2));
    H(1,2) = -sin(line_gl(2));
    H(1,3) = 0;
    % only non-zero at current line in global coordinate
    H(1,line_idx) = 1;
    H(1,line_idx+1) = robot_gl(1)*sin(line_gl(2))-robot_gl(2)*cos(line_gl(2));
    
    H(2,1) = 0;
    H(2,2) = 0;
    H(2,3) = -1;
    H(2,line_idx) = 0;
    H(2,line_idx+1) = 1;
else
    
    H(1,1) = cos(line_gl(2));
    H(1,2) = sin(line_gl(2));
    H(1,3) = 0;
    % only non-zero at current line in global coordinate
    H(1,line_idx) = -1;
    H(1,line_idx+1) = -robot_gl(1)*sin(line_gl(2))+robot_gl(2)*cos(line_gl(2));
    
    H(2,1) = 0;
    H(2,2) = 0;
    H(2,3) = -1;
    H(2,line_idx) = 0;
    H(2,line_idx+1) = 1;    
end


