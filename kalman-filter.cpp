#include <iostream>
#include "kalman-filter.h"

// void KalmanFilter::state_propagation(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_minus,\
//      Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_minus_k_minus, Eigen::MatrixXf R_xytheta, Eigen::VectorXf u, float dt) {
//          int state_num = state_all.rows();
//          Eigen::MatrixXf Fx(3,state_num);
//          Fx.setZero(3,state_num);
//          Fx(0,0) = 1;
//          Fx(1,1) = 1;
//          Fx(2,2) = 1;
// 
//          // 1.mean propagation
//          Eigen::Vector3f push_term;
//          push_term << robot_r/2*(u(1,0)+u(0,0))*cos(state_all(2,0))*dt,robot_r/2*(u(1,0)+u(0,0))*sin(state_all(2,0))*dt, robot_r/robot_b*(u(1,0)-u(0,0))*dt;
//          state_all_new = state_all + Fx.transpose()*push_term;
// 
//          // 2.covariance propagation
//          Eigen::MatrixXf J_push(3,3);
//          J_push.setZero(3,3);
//          J_push(0,2) = -robot_r/2*(u(1,0)+u(0,0))*sin(state_all(2,0))*dt;
//          J_push(1,2) = robot_r/2*(u(1,0)+u(0,0))*cos(state_all(2,0))*dt;
//          Eigen::MatrixXf G(state_num,state_num);
//          G = Eigen::MatrixXf::Identity(state_num,state_num)+Fx.transpose()*J_push*Fx;
//          Cov_k_k_minus = G*Cov_k_minus_k_minus*G.transpose()+Fx.transpose()*R_xytheta*Fx;
// }

void KalmanFilter::state_propagation(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_minus,\
     Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_minus_k_minus, Eigen::MatrixXf R_vw, Eigen::VectorXf u, float dt) {
           
         int state_num = state_all.rows();
         Eigen::MatrixXf Fx(3,state_num);
         Fx.setZero(3,state_num);
         Fx(0,0) = 1;
         Fx(1,1) = 1;
         Fx(2,2) = 1;

         // 1.mean propagation
         Eigen::Vector3f push_term;
         push_term << u(0,0)*cos(state_all(2,0))*dt,u(0,0)*sin(state_all(2,0))*dt, u(1,0)*dt;
         state_all_new = state_all + Fx.transpose()*push_term;

         // 2.covariance propagation
         Eigen::MatrixXf J_push(3,3);
         J_push.setZero(3,3);
         J_push(0,2) = -u(0,0)*sin(state_all(2,0))*dt;
         J_push(1,2) = u(0,0)*cos(state_all(2,0))*dt;
         Eigen::MatrixXf Phi(state_num,state_num);
         Phi = Eigen::MatrixXf::Identity(state_num,state_num)+Fx.transpose()*J_push*Fx;
         Eigen::MatrixXf G(3,2);
         G.setZero(3,2);
         G(0,0) = -dt*cos(state_all(2,0));
         G(1,0) = -dt*sin(state_all(2,0));
         G(2,1) = -dt;
         Cov_k_k_minus = Phi*Cov_k_minus_k_minus*Phi.transpose()+Fx.transpose()*G*R_vw*G.transpose()*Fx;
         //if(Cov_k_k_minus(0,0)<1.){
         //    Cov_k_k_minus(0,0) = 1.;
         //}
         //if(Cov_k_k_minus(1,1)<1.){
         //    Cov_k_k_minus(1,1)=1.;
         //}
         //float temp_cov = std::pow(0.00175,2);
         //if(Cov_k_k_minus(2,2)<temp_cov){
         //    Cov_k_k_minus(2,2)=temp_cov;
         //}
         

}

void KalmanFilter::measurement_update(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new,\
    Eigen::VectorXf state_all, Eigen::MatrixXf Cov_k_k_minus, Eigen::MatrixXf line_robot, Eigen::MatrixXf R_line, Eigen::VectorXi line_idx) {
        int num_state = state_all.rows();
        int num_match = line_idx.rows();


        Eigen::MatrixXf H_all(2*num_match,num_state);

        Eigen::VectorXf measure_dif_all(2*num_match);

        Eigen::MatrixXf R_all(2*num_match,2*num_match);
        R_all.setZero(2*num_match,2*num_match);
        for (int i = 0;i<num_match;i++){
            Eigen::Vector3f robot_gl;
            robot_gl << state_all(0,0),state_all(1,0),state_all(2,0);
            Eigen::Vector2f line_gl_single;
            line_gl_single << state_all(line_idx(i),0),state_all(line_idx(i)+1,0);
            int line_idx_single;
            line_idx_single = line_idx(i,0);

            //TODO: Not sure how to call this function
            int eq_idx = -1;
            Eigen::MatrixXf H_single;
            H_single.setZero(2,num_state);
            cal_H_mat(H_single, eq_idx,\
                        robot_gl, line_gl_single, line_idx_single, num_state);

//            std::cout << "H_single: \n" << H_single << std::endl;
//            std::cout << "eq_idx " << eq_idx;
            Eigen::Vector2f line_robot_single;
            line_robot_single = line_robot.col(i);
            Eigen::Vector2f measure_dif;
            if (eq_idx == 0){
                float orien_dif = line_robot_single(1,0)-line_gl_single(1,0)+robot_gl(2,0);
                measure_dif(0,0) = line_robot_single(0,0)-line_gl_single(0,0)+\
                                     robot_gl(0,0)*cos(line_gl_single(1,0))+robot_gl(1,0)*sin(line_gl_single(1,0));
                measure_dif(1,0) = atan2(sin(orien_dif),cos(orien_dif));

            }
            else{
                float orien_dif = line_robot_single(1,0)-line_gl_single(1,0)+robot_gl(2,0)-M_PI;
                measure_dif(0,0) = line_robot_single(0,0)+line_gl_single(0,0)\
                                     -robot_gl(0,0)*cos(line_gl_single(1,0))-robot_gl(1,0)*sin(line_gl_single(1,0));
                measure_dif(1,0) = atan2(sin(orien_dif),cos(orien_dif));
            }
            //TODO: Not sure how to add single matrix to the whole matrix and block ones

            H_all.row(2*i) = H_single.row(0);
            H_all.row(2*i+1) = H_single.row(1);
            measure_dif_all(2*i,0) = measure_dif(0,0);
            measure_dif_all(2*i+1,0) = measure_dif(1,0);
            R_all.block(2*i,2*i,2,2) = R_line.block(2*i,0,2,2);

        }
        // do the EKF measurement update
        Eigen::MatrixXf S(2*num_match,2*num_match);
        S = H_all*Cov_k_k_minus*H_all.transpose() + R_all;
        Eigen::MatrixXf K(num_state,2*num_match);
        K = Cov_k_k_minus*H_all.transpose()*S.inverse();
        // update the covariance
        Cov_k_k_new = Cov_k_k_minus-K*H_all*Cov_k_k_minus;
        // enforce symmetric
        Cov_k_k_new = (Cov_k_k_new+Cov_k_k_new.transpose())/2;
        // update the estimated variables
        state_all_new = state_all+K*measure_dif_all;

}

void KalmanFilter::update_mached_line_endpoint(Eigen::MatrixXf& endpoints_in_state_new,\
    Eigen::VectorXf state_all, Eigen::MatrixXf endpoints_in_state, Eigen::VectorXi line_idx, Eigen::MatrixXf endpoints_line_measure) {
        // 1. Extrac line parameter from state vector
        int mached_num = line_idx.rows();
        Eigen::MatrixXf line_state(2,mached_num);
        Eigen::MatrixXf endpoint_ref_1(2,mached_num);
        Eigen::MatrixXf endpoint_ref_2(2,mached_num);
        for (int i = 0; i<mached_num; i++){
            int idx = (line_idx(i)-3)/2;
            line_state.col(i) <<state_all(line_idx(i),0),state_all(line_idx(i)+1,0);
            endpoint_ref_1.col(i) << endpoints_in_state(0,idx), endpoints_in_state(2,idx);
            endpoint_ref_2.col(i) << endpoints_in_state(1,idx), endpoints_in_state(3,idx);
        }


        // 2. Change the measured lines ending points to global coordinate
        Eigen::Vector3f robot_vector;
        robot_vector << state_all(0,0),state_all(1,0),state_all(2,0);
        Eigen::MatrixXf end1_test_gl(2,mached_num);
        Eigen::MatrixXf end2_test_gl(2,mached_num);

        for (int i = 0; i<mached_num; i++){
            Eigen::Vector2f L_endpoints;
            L_endpoints << endpoints_line_measure(0,i), endpoints_line_measure(2,i);
            Eigen::Vector2f G_endpoints;
            G_endpoints = end_point_2_global(robot_vector, L_endpoints);
            end1_test_gl.col(i) = G_endpoints;
            L_endpoints << endpoints_line_measure(1,i), endpoints_line_measure(3,i);
            G_endpoints = end_point_2_global(robot_vector, L_endpoints);
            end2_test_gl.col(i) = G_endpoints;
        }


        // 3. Project all the ending points to the current lines due to the measurement update
        // 3.1 update the ending points in state
        int state_num = state_all.rows();
        Eigen::VectorXf whole_line_state(state_num-3);
        whole_line_state = state_all.middleRows(3, state_num-3);
        Eigen::MatrixXf endpoints_in_state_update;
        endpoints_in_state_update.setZero(4,endpoints_in_state.cols());
        p_2_line_state( endpoints_in_state_update,\
                        whole_line_state, endpoints_in_state);

//        std::cout << "endpoints_in_state_update:\n " << endpoints_in_state_update << std::endl;
//        std::cout << "whole_line_state:\n " << whole_line_state << std::endl;
//        std::cout << "endpoints_in_state:\n " << endpoints_in_state << std::endl;

        // 3.2 update the matched lines ending points in state and measurements
        Eigen::MatrixXf state_end_1;
        state_end_1.setZero(2,mached_num);
        Eigen::MatrixXf state_end_2;
        state_end_2.setZero(2,mached_num);
        Eigen::MatrixXf measure_end_1;
        measure_end_1.setZero(2,mached_num);
        Eigen::MatrixXf measure_end_2;
        measure_end_2.setZero(2,mached_num);
        p_2_line_match( state_end_1, state_end_2, measure_end_1, measure_end_2,\
                  line_state, endpoint_ref_1, endpoint_ref_2, end1_test_gl, end2_test_gl);


        // 4.Update the ending points
        endpoints_in_state_new = endpoints_in_state_update;
        for (int i = 0; i<mached_num; i++){
            // get the line index that needs to update the end points
            int idx = (line_idx(i)-3)/2;
            float theta_ref = line_state(1,i);
            // compute the relative angle w.r.t. the line parameter theta
            float end1_test_theta = atan2(measure_end_1(1,i),measure_end_1(0,i))-theta_ref;
            end1_test_theta = atan2(sin(end1_test_theta),cos(end1_test_theta));

            float end2_test_theta = atan2(measure_end_2(1,i),measure_end_2(0,i))-theta_ref;
            end2_test_theta = atan2(sin(end2_test_theta),cos(end2_test_theta));

            float end1_ref_theta = atan2(state_end_1(1,i),state_end_1(0,i))-theta_ref;
            end1_ref_theta = atan2(sin(end1_ref_theta),cos(end1_ref_theta));

            float end2_ref_theta = atan2(state_end_2(1,i),state_end_2(0,i))-theta_ref;
            end2_ref_theta = atan2(sin(end2_ref_theta),cos(end2_ref_theta));

            if(end1_test_theta > end2_test_theta){
                float temp = end1_test_theta;
                end1_test_theta = end2_test_theta;
                end2_test_theta = temp;
                Eigen::Vector2f temp_vec;
                temp_vec = measure_end_1.col(i);
                measure_end_1.col(i) = measure_end_2.col(i);
                measure_end_2.col(i) = temp_vec;

            }
            if(end1_ref_theta > end2_ref_theta){
                float temp = end1_ref_theta;
                end1_ref_theta = end2_ref_theta;
                end2_ref_theta = temp;
                Eigen::Vector2f temp_vec;
                temp_vec = state_end_1.col(i);
                state_end_1.col(i) = state_end_2.col(i);
                state_end_2.col(i) = temp_vec;
            }
            if((end1_test_theta >= end1_ref_theta)&&(end2_test_theta<=end2_ref_theta)){
            }
            else if((end1_test_theta>=end1_ref_theta)&&(end1_test_theta<=end2_ref_theta)){
                state_end_2.col(i) = measure_end_2.col(i);
            }
            else if((end2_test_theta>=end1_ref_theta)&&(end2_test_theta<=end2_ref_theta)){
                state_end_1.col(i) = measure_end_1.col(i);
            }
            else if((end1_ref_theta>= end1_test_theta)&&(end2_ref_theta<= end2_test_theta)&&(end2_ref_theta >= end1_test_theta)){
                state_end_1.col(i) = measure_end_1.col(i);
                state_end_2.col(i) = measure_end_2.col(i);
            }
            endpoints_in_state_new.col(idx)<< state_end_1(0,i),state_end_2(0,i),state_end_1(1,i),state_end_2(1,i);

        }

}

void KalmanFilter::cov_update_new_line_measure(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new,\
    Eigen::MatrixXf Cov_k_k, Eigen::Matrix2f R_new_line, Eigen::VectorXf state_all, Eigen::Vector2f new_line_para) {
        // 1.determine which measurement equation should be used
        float d_line = new_line_para(0,0);
        float theta_line = new_line_para(1,0);
        Eigen::Vector2f d_gl_candi;
        d_gl_candi(0,0) = d_line+state_all(0,0)*cos(theta_line+state_all(2,0))+state_all(1,0)*sin(theta_line+state_all(2,0));
        d_gl_candi(1,0) = -d_line+state_all(0,0)*cos(theta_line+state_all(2,0)-M_PI)+state_all(1,0)*sin(theta_line+state_all(2,0)-M_PI);
        int eq_use,c;
        float s = d_gl_candi.maxCoeff(&eq_use, &c);
        Eigen::Vector2f theta_gl_candi;
        theta_gl_candi << theta_line+state_all(2,0),theta_line+state_all(2,0)-M_PI;

        // 2.construct the new line state
        Eigen::Vector2f new_line_state;
        new_line_state << d_gl_candi(eq_use,0),atan2(sin(theta_gl_candi(eq_use,0)),cos(theta_gl_candi(eq_use,0)));

        // 3.calculate the Jacobian matrix
        Eigen::MatrixXf J_state;
        J_state.setZero(2,state_all.rows());
        Eigen::Matrix2f J_line;
        J_line.setZero(2,2);
        if(eq_use==0){
            // Jacobian w.r.t. the state vector
            J_state(0,0) = cos(theta_line+state_all(2,0));
            J_state(0,1) = sin(theta_line+state_all(2,0));
            J_state(0,2) = -state_all(0,0)*sin(theta_line+state_all(2,0))+state_all(1,0)*cos(theta_line+state_all(2,0));
            J_state(1,0) = 0;
            J_state(1,1) = 0;
            J_state(1,2) = 1;

            // Jacobian w.r.t. the new line parameter
            J_line(0,0) = 1;
            J_line(0,1) = -state_all(0,0)*sin(theta_line+state_all(2,0))+state_all(1,0)*cos(theta_line+state_all(2,0));
            J_line(1,0) = 0; 
            J_line(1,1) = 1;
        }
        else{
            // Jacobian w.r.t. the state vector
            J_state(0,0) = cos(theta_line+state_all(2,0)-M_PI);
            J_state(0,1) = sin(theta_line+state_all(2,0)-M_PI);
            J_state(0,2) = -state_all(0,0)*sin(theta_line+state_all(2,0)-M_PI)+state_all(1,0)*cos(theta_line+state_all(2,0)-M_PI);
            J_state(1,0) = 0;
            J_state(1,1) = 0;
            J_state(1,2) = 1;

            // Jacobian w.r.t. the new line parameter
            J_line(0,0) = -1;
            J_line(0,1) = -state_all(0,0)*sin(theta_line+state_all(2,0)-M_PI)+state_all(1,0)*cos(theta_line+state_all(2,0)-M_PI);
            J_line(1,0) = 0;
            J_line(1,1) = 1; 
        }

        // 3. update the state to include new line
        state_all_new.setZero(state_all.rows()+2,1);
        state_all_new.topRows(state_all.rows()) = state_all;
        state_all_new.bottomRows(2) = new_line_state;

        // 4. update the state covariance
        Cov_k_k_new.setZero(state_all.rows()+2,state_all.rows()+2);
        Cov_k_k_new.block(0, 0, state_all.rows(), state_all.rows()) = Cov_k_k;
        Cov_k_k_new.block(0,state_all.rows(),state_all.rows(),2) = Cov_k_k*J_state.transpose();
        Cov_k_k_new.block(state_all.rows(), 0, 2, state_all.rows()) = J_state*Cov_k_k;
        Cov_k_k_new.block(state_all.rows(), state_all.rows(), 2, 2) = J_line*R_new_line*J_line.transpose();

    }

void KalmanFilter::cal_H_mat(Eigen::MatrixXf& H, int& eq_idx,\
    Eigen::Vector3f robot_gl,Eigen::Vector2f line_gl, int line_idx, int num_state) {

        //decide which measurement equation should be used
        Eigen::Vector2f d_line;
        d_line<<line_gl(0,0)-robot_gl(0,0)*cos(line_gl(1,0))-robot_gl(1,0)*sin(line_gl(1,0)),\
        -line_gl(0,0)+robot_gl(0,0)*cos(line_gl(1,0))+robot_gl(1,0)*sin(line_gl(1,0));
        int c;
        float s = d_line.maxCoeff(&eq_idx, &c);
        H.setZero(2,num_state);

        if(eq_idx==0){
            H(0,0) = -cos(line_gl(1,0));
            H(0,1) = -sin(line_gl(1,0));
            H(0,2) = 0;
            //only non-zero at current line in global coordinate
            H(0,line_idx) = 1;
            H(0,line_idx+1) = robot_gl(0,0)*sin(line_gl(1,0))-robot_gl(1,0)*cos(line_gl(1,0));

            H(1,0) = 0;
            H(1,1) = 0;
            H(1,2) = -1;
            H(1,line_idx) = 0;
            H(1,line_idx+1) = 1;

        }
        else {
            H(0, 0) = cos(line_gl(1, 0));
            H(0, 1) = sin(line_gl(1, 0));
            H(0, 2) = 0;
            //% only non-zero at current line in global coordinate
            H(0, line_idx) = -1;
            H(0, line_idx + 1) = -robot_gl(0,0) * sin(line_gl(1,0)) + robot_gl(1,0) * cos(line_gl(1,0));

            H(1, 0) = 0;
            H(1, 1) = 0;
            H(1, 2) = -1;
            H(1, line_idx) = 0;
            H(1, line_idx + 1) = 1;

        }
}

void KalmanFilter::add_line_2_state(Eigen::VectorXf& state_all_new, Eigen::MatrixXf& Cov_k_k_new, Eigen::Vector4f& end_point_gl,\
     Eigen::VectorXf state_all, Eigen::MatrixXf Cov, Eigen::Vector2f x_part, Eigen::Vector2f y_part, Eigen::Matrix2f R_add_single) {

         Eigen::Vector3f robot_vector;
         robot_vector<< state_all(0,0),state_all(1,0),state_all(2,0);
         Eigen::Vector2f L_endpoints_1;
         L_endpoints_1 << x_part(0,0), y_part(0,0);
         Eigen::Vector2f G_endpoints_1;
         G_endpoints_1 = end_point_2_global(robot_vector, L_endpoints_1);
         Eigen::Vector2f L_endpoints_2;
         L_endpoints_2 << x_part(1,0), y_part(1,0);
         Eigen::Vector2f G_endpoints_2;
         G_endpoints_2 = end_point_2_global(robot_vector, L_endpoints_2);
         end_point_gl<<G_endpoints_1(0,0),G_endpoints_2(0,0),G_endpoints_1(1,0),G_endpoints_2(1,0);

         Eigen::Vector2f end_point_1_r;
         end_point_1_r = end_point_2_robot(L_endpoints_1);
         Eigen::Vector2f end_point_2_r;
         end_point_2_r = end_point_2_robot(L_endpoints_2);  

         Eigen::Vector2f d_est;
         d_est << end_point_1_r.norm(), end_point_2_r.norm();
         Eigen::Vector2f theta_est;
         theta_est << atan2(end_point_1_r(1,0),end_point_1_r(0,0)),atan2(end_point_2_r(1,0),end_point_2_r(0,0));
         float d_l,theta_l;
         line_para_cal(d_est, theta_est, d_l, theta_l);

         Eigen::Vector2f line_add_single_r;
         line_add_single_r << d_l, theta_l;
         //TODO cov_update function

         cov_update_new_line_measure(state_all_new, Cov_k_k_new,\
          Cov, R_add_single, state_all, line_add_single_r);




}
