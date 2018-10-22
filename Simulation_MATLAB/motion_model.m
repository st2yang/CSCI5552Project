function x_new = motion_model(state_all,u,dt)


%push_term = [r/2*(u(2)+u(1))*cos(state_all(3))*dt,r/2*(u(2)+u(1))*sin(state_all(3))*dt, r/b*(u(2)-u(1))*dt]';
push_term = [u(1)*dt*cos(state_all(3)), u(1)*dt*sin(state_all(3)), u(2)*dt]';

x_new = state_all+push_term;