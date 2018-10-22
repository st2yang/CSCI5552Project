
function [state_all_new,Cov_k_k_new] = cov_update_new_line_measure(Cov_k_k,R_new_line,state_all,new_line_para)

% This function is used to calculate the covariance matrix with newly-added
% lines as another new states. 
% This function should be performed after all the measurement updates are
% done. The remaining new lines will be added to the state, as a result, we
% need to calculate the covariance matrix again.

% Input:
%       Cov_k_k: (3+2n)x(3+2n) covariance matrix k|k
%       R_new_line: 2x2 covariance matrix of the estimated new line 
%       state_all: (3+2n)x1 vector, robot x,y,theta are at the first 3
%       new_line_para: 2x1 vector, new line's parameter,(d,theta)',this is
%       in robot coordinate


% 1.determine which measurement equation should be used
d_line = new_line_para(1);
theta_line = new_line_para(2);
d_gl_candi = zeros(2,1);
d_gl_candi(1) = d_line+state_all(1)*cos(theta_line+state_all(3))...
    +state_all(2)*sin(theta_line+state_all(3));
d_gl_candi(2) = -d_line+state_all(1)*cos(theta_line+state_all(3)-pi)...
    +state_all(2)*sin(theta_line+state_all(3)-pi);
[~,eq_use] = max(d_gl_candi);
theta_gl_candi = zeros(2,1);
theta_gl_candi(1) = theta_line+state_all(3);
theta_gl_candi(2) = theta_line+state_all(3)-pi;

% 2.construct the new line state
new_line_state = zeros(2,1);
new_line_state(1) = d_gl_candi(eq_use);
new_line_state(2) = atan2(sin(theta_gl_candi(eq_use)),cos(theta_gl_candi(eq_use)));

% 3.calculate the Jacobian matrix
J_state = zeros(2,length(state_all));
J_line = zeros(2);
if eq_use==1
    % Jacobian w.r.t. the state vector
    J_state(1,1) = cos(theta_line+state_all(3));
    J_state(1,2) = sin(theta_line+state_all(3));
    J_state(1,3) = -state_all(1)*sin(theta_line+state_all(3))+...
        state_all(2)*cos(theta_line+state_all(3));
    J_state(2,1) = 0;
    J_state(2,2) = 0;
    J_state(2,3) = 1;
    
    % Jacobian w.r.t. the new line parameter
    J_line(1,1) = 1;
    J_line(1,2) = -state_all(1)*sin(theta_line+state_all(3))+...
        state_all(2)*cos(theta_line+state_all(3));
    J_line(2,1) = 0;
    J_line(2,2) = 1;
else
    % Jacobian w.r.t. the state vector
    J_state(1,1) = cos(theta_line+state_all(3)-pi);
    J_state(1,2) = sin(theta_line+state_all(3)-pi);
    J_state(1,3) = -state_all(1)*sin(theta_line+state_all(3)-pi)+...
        state_all(2)*cos(theta_line+state_all(3)-pi);
    J_state(2,1) = 0;
    J_state(2,2) = 0;
    J_state(2,3) = 1;
    
    % Jacobian w.r.t. the new line parameter
    J_line(1,1) = -1;
    J_line(1,2) = -state_all(1)*sin(theta_line+state_all(3)-pi)+...
        state_all(2)*cos(theta_line+state_all(3)-pi);
    J_line(2,1) = 0;
    J_line(2,2) = 1; 
end

% 3. update the state to include new line
state_all_new = [state_all;new_line_state];
% 4. update the state covariance
Cov_k_k_new = [Cov_k_k, Cov_k_k*J_state';J_state*Cov_k_k, J_line*R_new_line*J_line'];

end
    
    
    