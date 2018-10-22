
function [state_all_new,Cov_k_k_minus] = state_propagation(state_all,Cov_k_minus_k_minus,R_vw,u,dt)

% This function do the state propagation based on the motion model
% considered in the paper "modelling of mobile robot dynamics".

% Input:
%       state_all: (3+2n)x1 vector,[x,y,theta,d1,theta1,...]'
%       Cov_k_minus_k_minus: (3+2n)x(3+2n) covariance matrix
%       R_xytheta: 3x3 covariance matrix for robot's x,y,theta
%       u: 2x1 input vector,[left_angular_velocity,
%                            right_angular_velocity]'
%       dt: sampling time
%       r: robot wheel parameter
%       b: robot length parameter

Fx = zeros(3,length(state_all));
Fx(1:3,1:3) = eye(3);

% 1.mean propagation
%push_term = [r/2*(u(2)+u(1))*cos(state_all(3))*dt,r/2*(u(2)+u(1))*sin(state_all(3))*dt, r/b*(u(2)-u(1))*dt]';
push_term = [u(1)*dt*cos(state_all(3)), u(1)*dt*sin(state_all(3)), u(2)*dt]';
state_all_new = state_all + Fx'*push_term;


% 2.covariance propagation
J_push = zeros(3);
%J_push(1,3) = -r/2*(u(2)+u(1))*sin(state_all(3))*dt;
%J_push(2,3) = r/2*(u(2)+u(1))*cos(state_all(3))*dt;
J_push(1,3) = -u(1)*dt*sin(state_all(3));
J_push(2,3) = u(1)*dt*cos(state_all(3));
Phi = eye(length(state_all))+Fx'*J_push*Fx;
G = [-dt*cos(state_all(3)),0;-dt*sin(state_all(3)) 0; 0 -dt];

Cov_k_k_minus = Phi*Cov_k_minus_k_minus*Phi'+Fx'*G*R_vw*G'*Fx;

end
