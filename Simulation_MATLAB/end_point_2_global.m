function    [end_point_1,end_point_2,end_point_gl] = end_point_2_global(robot_vector,x_part,y_part)


theta_robot = robot_vector(3);
x_robot = robot_vector(1);
y_robot = robot_vector(2);
Rot_rob_2_gl = [cos(theta_robot), -1*sin(theta_robot), x_robot;sin(theta_robot), cos(theta_robot), y_robot;0 0 1];
Rot_line_2_rob = [cos(-pi/2), -1*sin(-pi/2) 0;sin(-pi/2), cos(-pi/2), 0;0 0 1];
end_point_1_temp = Rot_rob_2_gl*Rot_line_2_rob*[x_part(1), y_part(1), 1]';
end_point_1 = end_point_1_temp(1:2,1);
end_point_2_temp = Rot_rob_2_gl*Rot_line_2_rob*[x_part(2), y_part(2), 1]';
end_point_2 = end_point_2_temp(1:2,1);
end_point_gl = [end_point_1(1),end_point_2(1),end_point_1(2),end_point_2(2)]';

