function plot_endpoints_state(endpoints_in_state)

for i = 1:size(endpoints_in_state,2)
    end_point_gl = endpoints_in_state(:,i);
    plot([end_point_gl(1),end_point_gl(2)]'/1000,[end_point_gl(3),end_point_gl(4)]'/1000,'r-','LineWidth',8);  
end