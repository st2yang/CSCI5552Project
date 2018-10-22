#include "line-match.h"

void LineMatch::performLineMatching(Eigen::MatrixXf& line_robot, Eigen::MatrixXf& R_line, Eigen::VectorXi& line_idx, Eigen::MatrixXf& endpoints_line_measure, Eigen::MatrixXf& line_un_match, Eigen::MatrixXf& endpoints_un_match, Eigen::MatrixXf& R_un_match,\
  Eigen::MatrixXf line_in_state, Eigen::MatrixXf endpoints_in_state, Eigen::MatrixXf line_test, Eigen::MatrixXf endpoints_test ,Eigen::MatrixXf R_test, Eigen::Vector3f robot_vector) {

    Eigen::Vector2f trans_move;
    trans_move << robot_vector(0,0),robot_vector(1,0);
    float theta_turn = robot_vector(2,0);

    const int num_line_state = line_in_state.cols();
    const int num_line_test = line_test.cols();

    int match_flag = 0;

    Eigen::Vector2f G_endpoints;
    Eigen::Vector2f L_endpoints;

    Eigen::MatrixXf end1_test_gl(2,num_line_test);
    Eigen::MatrixXf end2_test_gl(2,num_line_test);

    Eigen::Vector2f R_endpoint_1;
    Eigen::Vector2f R_endpoint_2;

    Eigen::Vector2f d_est;
    Eigen::Vector2f theta_est;

    Eigen::MatrixXf line_test_r(2,num_line_test);


    for(int k = 0; k < num_line_test; k++) {
        L_endpoints << endpoints_test(0,k), endpoints_test(2,k);
        G_endpoints = end_point_2_global(robot_vector,L_endpoints);
        end1_test_gl.col(k) << G_endpoints;
        R_endpoint_1 = end_point_2_robot(L_endpoints);

        L_endpoints << endpoints_test(1,k), endpoints_test(3,k);
        G_endpoints = end_point_2_global(robot_vector,L_endpoints);
        end2_test_gl.col(k) << G_endpoints;
        R_endpoint_2 = end_point_2_robot(L_endpoints);
        d_est << sqrt(R_endpoint_1.transpose()*R_endpoint_1),sqrt(R_endpoint_2.transpose()*R_endpoint_2);
        theta_est<< atan2(R_endpoint_1(1),R_endpoint_1(0)),atan2(R_endpoint_2(1),R_endpoint_2(0));

        float d_l;
        float theta_l;
        line_para_cal(d_est,theta_est, d_l, theta_l);
        line_test_r.col(k)<<d_l,theta_l;
    }

    Eigen::MatrixXf pair_val(num_line_state,num_line_test);

    for(int k = 0;k<num_line_state;k++){
        Eigen::Vector2f line_single;
        line_single = line_in_state.col(k);
        float theta_ref;
        theta_ref = line_single(1,0);
        float d_ref;
        d_ref = line_single(0,0);
        Eigen::Vector4f endpoint_ref;
        endpoint_ref = endpoints_in_state.col(k);
        for(int j=0;j<num_line_test;j++){
            Eigen::Matrix<float,2,2> R_single;
            R_single = R_test.block(2*j,0,2,2);
            float theta_test;
            theta_test = line_test_r(1,j);

            // orientation difference
            Eigen::Vector2f dif_orien;
            dif_orien<<theta_test-theta_ref+theta_turn-M_PI,theta_test-theta_ref+theta_turn;
            Eigen::Vector2f d_test_m;
            d_test_m<< -d_ref+trans_move(0,0)*cos(theta_ref)+trans_move(1,0)*sin(theta_ref),\
            d_ref-trans_move(0,0)*cos(theta_ref)-trans_move(1,0)*sin(theta_ref);
            int measure_eq_idx,c;
            float s = d_test_m.maxCoeff(&measure_eq_idx, &c);
            dif_orien << abs(atan2(sin(dif_orien(0,0)), cos(dif_orien(0,0)))),\
             abs(atan2(sin(dif_orien(1,0)), cos(dif_orien(1,0))));
            float orien_com_val;
            orien_com_val = dif_orien(measure_eq_idx,0);

            // compute the test line midpoint distance to the reference line
            Eigen::Vector2f end1_test_gl_single;
            end1_test_gl_single = end1_test_gl.col(j);
            Eigen::Vector2f end2_test_gl_single;
            end2_test_gl_single = end2_test_gl.col(j);
            Eigen::Vector2f mid_test_gl;
            mid_test_gl = end2_test_gl_single+0.5*(end1_test_gl_single-end2_test_gl_single);
            float dist_com_val;
            dist_com_val = abs(cos(theta_ref)*mid_test_gl(0,0)+sin(theta_ref)*mid_test_gl(1,0)-d_ref);

            // compute the overlapping rate
            Eigen::Vector2f end1_test_gl_in_line;
            end1_test_gl_in_line(0,0) = sin(theta_ref)*(sin(theta_ref)*end1_test_gl_single(0,0)-cos(theta_ref)*end1_test_gl_single(1,0))\
            -cos(theta_ref)*(-1*d_ref);
            end1_test_gl_in_line(1,0) = cos(theta_ref)*(-1*sin(theta_ref)*end1_test_gl_single(0,0)+cos(theta_ref)*end1_test_gl_single(1,0))\
            -sin(theta_ref)*(-1*d_ref);
            float end1_test_theta;
            end1_test_theta = atan2(end1_test_gl_in_line(1,0),end1_test_gl_in_line(0,0))-theta_ref;
            end1_test_theta = atan2(sin(end1_test_theta),cos(end1_test_theta));

            Eigen::Vector2f end2_test_gl_in_line;
            end2_test_gl_in_line(0,0) = sin(theta_ref)*(sin(theta_ref)*end2_test_gl_single(0,0)-cos(theta_ref)*end2_test_gl_single(1,0))\
            -cos(theta_ref)*(-1*d_ref);
            end2_test_gl_in_line(1,0) = cos(theta_ref)*(-1*sin(theta_ref)*end2_test_gl_single(0,0)+cos(theta_ref)*end2_test_gl_single(1,0))\
            -sin(theta_ref)*(-1*d_ref);
            float end2_test_theta;
            end2_test_theta = atan2(end2_test_gl_in_line(1,0),end2_test_gl_in_line(0,0))-theta_ref;
            end2_test_theta = atan2(sin(end2_test_theta),cos(end2_test_theta));

            float end1_ref_theta;
            end1_ref_theta = atan2(endpoint_ref(2,0),endpoint_ref(0,0))-theta_ref;
            end1_ref_theta = atan2(sin(end1_ref_theta),cos(end1_ref_theta));

            float end2_ref_theta;
            end2_ref_theta = atan2(endpoint_ref(3,0),endpoint_ref(1,0))-theta_ref;
            end2_ref_theta = atan2(sin(end2_ref_theta),cos(end2_ref_theta));

            if (end1_test_theta > end2_test_theta){
                float temp;
                temp = end1_test_theta;
                end1_test_theta = end2_test_theta;
                end2_test_theta = temp;
            }
            if (end1_ref_theta > end2_ref_theta){
                float temp;
                temp = end1_ref_theta;
                end1_ref_theta = end2_ref_theta;
                end2_ref_theta = temp;

            }
            float lap_rate;
            if((end1_test_theta >= end1_ref_theta)&& (end2_test_theta<=end2_ref_theta)){
                lap_rate = 1;
            }
            else if((end1_test_theta>=end1_ref_theta)&&(end1_test_theta<=end2_ref_theta)){
                lap_rate = (end2_ref_theta-end1_test_theta)/(end2_ref_theta-end1_ref_theta);

            }
            else if((end2_test_theta>=end1_ref_theta)&&(end2_test_theta<=end2_ref_theta)){
                lap_rate = (end2_test_theta-end1_ref_theta)/(end2_ref_theta-end1_ref_theta);
            }
            else if((end1_ref_theta>= end1_test_theta)&&(end2_ref_theta<= end2_test_theta)&&(end2_ref_theta >= end1_test_theta)){
                lap_rate = 1;
            }
            else{
                lap_rate = 0;
            }

            if((orien_com_val<=orien_thresh)&&(dist_com_val<=dist_thresh)&&(lap_rate>=0.1)){
                pair_val(k,j) = 0.2*orien_com_val/orien_thresh+0.3*dist_com_val/dist_thresh+0.5*(1-lap_rate)/0.9;
                match_flag = 1;
            }
            else{
                pair_val(k,j) = 10;
            }

        }
    }

    // obtain the matching based on the pair_val
    if (match_flag>0){

        Eigen::MatrixXi pair_match(0,2);
        int count_match = 0;

        for(int i = 0;i<num_line_state;i++){
            //num_line_test
            Eigen::VectorXf pair_val_single(num_line_test);
            pair_val_single = pair_val.row(i);
            if(pair_val_single.minCoeff()<=1){
                int idx_match,c;
                float s = pair_val_single.minCoeff(&idx_match, &c);
                Eigen::MatrixXi temp(1,2);
                temp << i, idx_match;
                AddMatrixToRow(pair_match,temp);
                count_match++;
            }

        }

        //update the ending points for the matched lines in state
        int pair_match_state_num = pair_match.rows();
        Eigen::VectorXi temp(pair_match_state_num);
        line_robot.setZero(2,pair_match_state_num);
        R_line.setZero(2*pair_match_state_num,2);
        endpoints_line_measure.setZero(4,pair_match_state_num);
        line_idx=temp;
        for(int i = 0;i<pair_match_state_num;i++){
            int ref_idx = pair_match(i,0);
            int test_idx = pair_match(i,1);
            line_idx(i,0) = ref_idx*2+3;
            line_robot.col(i) = line_test_r.col(test_idx);
            R_line.block(2*i,0,2,2) = R_test.block(2*test_idx,0,2,2);
            endpoints_line_measure.col(i) = endpoints_test.col(test_idx);

        }

        line_un_match = line_test;
        endpoints_un_match = endpoints_test;
        R_un_match = R_test;

        Eigen::VectorXi idx_match_c2 = pair_match.col(1);
        Eigen::VectorXi idx_match = unique(idx_match_c2);

        for(int i = 0; i< idx_match.rows(); i++) {
            removeColumn(line_un_match,idx_match(i));
            removeColumn(endpoints_un_match,idx_match(i));
            removeRow(R_un_match,2*idx_match(i)+1);
            removeRow(R_un_match,2*idx_match(i));
        }

    }
 else {
        printf ( "Warning! No matching lines!" );
        exit ( 1 );
    }


}