#include <iostream>
#include <fstream>
#include <chrono>
#include "Aria.h"
#include <time.h>
#include <string>

#include "commons.h"
#include "kalman-filter.h"
#include "LineExtraction_RANSAC.h"
#include "line-match.h"
#include "movementcontroller.h"

//int main() {
//
//    string outfileName = "./data/endpoints.txt";
//    ofstream outFile(outfileName);
//
//    outFile << 1 << " " << 2;
//}


int main(int argc, char **argv) {

    // setup outfiles
    // odometry
    std::string endfileName = "/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints.txt";
    std::ofstream endFile(endfileName);

    std::string posefileName = "/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotpose.txt";
    std::ofstream poseFile(posefileName);

    std::string covfileName = "/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/robotcov.txt";
    std::ofstream covFile(covfileName);   
    
    std::string end_update_name = "/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/endpoints_update.txt";
    std::ofstream end_update_file(end_update_name); 

    std::string scan_name = "/home/marsyang/Documents/5552Folder/SLAM_DUNK/data/scan_update.txt";
    std::ofstream laser_file(scan_name); 
    

////////////////////////// Initialize Aria
    Aria::init();
    ArArgumentParser argParser(&argc, argv);
    argParser.loadDefaultArguments();
    ArRobot robot;
    ArSick sick;


    ArSimpleConnector connector(&argc, argv);

    if (!connector.parseArgs() || argc > 1)
    {
        connector.logOptions();
        exit(1);
    }

    ArKeyHandler keyHandler;
    Aria::setKeyHandler(&keyHandler);
    robot.attachKeyHandler(&keyHandler);

    ArSonarDevice sonar;
    robot.addRangeDevice(&sonar);
    robot.addRangeDevice(&sick);

    if (!connector.connectRobot(&robot))
    {
        printf("Could not connect to robot... exiting\n");
        Aria::shutdown();
        return 1;
    }

    robot.runAsync(true);


    //Initialize and connect to the laser
    sick.configureShort(false,ArSick::BAUD38400,ArSick::DEGREES180,ArSick::INCREMENT_ONE);
    connector.setupLaser(&sick);
    sick.runAsync();


    if (!sick.blockingConnect())
    {
        printf("Could not connect to SICK laser... exiting\n");
        robot.disconnect();
        Aria::shutdown();
        return 1;
    }


    //Enable the moters
    robot.enableMotors();
    robot.comInt(ArCommands::SOUNDTOG, 0);

    //Setup timers to manage dt loop timings
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    MovementController* mov = new MovementController(&robot, &sick);
    mov->start();

    robot.requestEncoderPackets();

////////////////////////// main part
    //Initialize origin
    Eigen::Vector3f robot_initial;
    robot_initial << 5/2.0*1000,0.0, M_PI/2;
    Eigen::Matrix3f Cov_initial;
    Cov_initial.setZero(3,3);
     Cov_initial(0,0) = std::pow(1,2);
     Cov_initial(1,1) = std::pow(1,2);
     Cov_initial(2,2) = std::pow(0.00175,2);

    Eigen::VectorXf state_all = robot_initial;
    Eigen::MatrixXf Cov = Cov_initial;

    // decalre tools here
    KalmanFilter* ekf = new KalmanFilter();
    LineMatch* line_match = new LineMatch();
    const std::list<ArSensorReading *> *readings;
    std::list<ArSensorReading *>::const_iterator it;
    Eigen::MatrixXf endpoints_in_state(0,0);

    // declare common variables here
    Eigen::MatrixXf line_un_match(0,0);
    Eigen::MatrixXf endpoints_un_match(0,0);
    Eigen::MatrixXf R_un_match(0,0);
    Eigen::MatrixXf line_add(0,0);
    Eigen::MatrixXf R_add(0,0);
    Eigen::MatrixXf R_add_single(0,0);
    Eigen::MatrixXf endpoints_add(0,0);
    Eigen::Vector2f x_part;
    Eigen::Vector2f y_part;
    Eigen::VectorXf state_all_new(1);
    Eigen::MatrixXf Cov_k_k_new(1,1);
    Eigen::Vector4f end_point_gl;
    float dt = 0.001;

    // declare const paras here
//     float sigma_u_x = 50;
//     float sigma_u_y = 50;
//     float sigma_u_theta = 5/180*3.14;
// 
//     Eigen::Matrix3f R_xytheta;
//     R_xytheta.setZero(3,3);
    
    float sigma_v = 10;
    float sigma_w = 2.5/180*3.14;
    Eigen::Matrix2f R_vw;
    R_vw.setZero(2,2);

//     Eigen::VectorXf encoder1(2,1);
//     encoder1 << 0,0;
//     Eigen::VectorXf encoder2(2,1);
//     encoder2 = encoder1;
    Eigen::VectorXf u(2,1);

    Eigen::VectorXf state_all_temp(1);
    Eigen::MatrixXf Cov_k_k_minus(1,1);



    Eigen::Matrix2f Cov_Laser;
    float sigma_d = 0.005;
    Cov_Laser << std::pow(2*sigma_d*1000,2),0,
            0,2*1e-8;

    float thresh = 250; //maximum distance from point to the line
    int min_num = 10;
    int maxIter = 1000;
    int min_rest = 20;
    float theta_thresh = 1;
    float seg_dist_thresh = 500;
    float line_seg_dist_thresh = 500;

    ////////////////////////////////////////////
    
    int k = 0;
    float loopTime = 0.0;
    while(1) {
            
      t1 = t2;
      t2 = std::chrono::high_resolution_clock::now();
      long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      dt = (float)microseconds / 1000000.0;
      

//         R_xytheta(0,0) = std::pow(2*dt*sigma_u_x,2)/10;
//         R_xytheta(1,1) = std::pow(2*dt*sigma_u_y,2)/10;
//         R_xytheta(2,2) = std::pow(2*dt*sigma_u_theta,2)/10;


//         encoder1 = encoder2;
//         encoder2(0) = robot.getLeftEncoder();
//         encoder2(1) = robot.getRightEncoder();
// 
//         u = (encoder2-encoder1)/65536*2*M_PI;

        robot.lock();
        float V = robot.getVel();
        float RTV = robot.getRotVel()* M_PI / 180.0;
        robot.unlock();
        u(0) = V;
        u(1) = RTV;

        R_vw(0,0) = std::pow(0.3*sigma_v,2);
        R_vw(1,1) = std::pow(0.5*sigma_w,2);


        //////// ekf propagate

        ekf->state_propagation(state_all_temp, Cov_k_k_minus,\
      state_all, Cov, R_vw, u, dt);
        std::cout << "Propagate_Cov:\n " << Cov_k_k_minus(0,0) << std::endl;
        state_all = state_all_temp;
        state_all(2,0) = atan2(sin(state_all(2,0)),cos(state_all(2,0)));
        Cov = Cov_k_k_minus;

        //////// laser and line extraction
        Eigen::MatrixXf d_est(0,0);
        Eigen::MatrixXf theta_est(0,0);
        Eigen::MatrixXf angle(1,1);
        Eigen::MatrixXf dist(1,1);

        sick.lockDevice();
        readings = sick.getRawReadings();
        sick.unlockDevice();

        for (it = readings->begin();it != readings->end();it++) {
            angle(0,0) = ((*it)->getSensorTh()+90.0)*M_PI/180.0;
            dist(0,0) = ((*it)->getRange());
//             std::cout << "angle: " << angle(0,0) << std::endl;
// 	    std::cout << "dist: " << dist(0,0) << std::endl;
            if (dist(0,0) > MAX_DIST) continue;
            AddMatrixToRow(d_est,dist);
            AddMatrixToRow(theta_est,angle);
        }
                

        LineExtraction_RANSAC line_extract;
        line_extract.Init(d_est,theta_est,Cov_Laser,thresh,min_num,maxIter,min_rest,theta_thresh, seg_dist_thresh, line_seg_dist_thresh);
        line_extract.LineExtraction_main();

        Eigen::MatrixXf line_para = line_extract.line;
        Eigen::MatrixXf endpoints_in_line = line_extract.endpoints_in_line;
        Eigen::MatrixXf R_mat = line_extract.R_mat;
       
        std::cout << "---------------------- \n";
        std::cout << "line_para:\n " << line_para << std::endl;
        std::cout << "endpoints_in_line:\n " << endpoints_in_line << std::endl;

        //////// line matching
        if (state_all.size()>3) {
            Eigen::VectorXf state_line = state_all.bottomRows(state_all.size()-3);
            Eigen::MatrixXf line_in_state = state_line;
            line_in_state.resize(2,state_line.size()/2);
            Eigen::MatrixXf endpoints_in_state_use = endpoints_in_state;
            Eigen::MatrixXf line_test = line_para;
            Eigen::MatrixXf endpoints_test = endpoints_in_line;
            Eigen::MatrixXf R_test = R_mat;
            Eigen::Vector3f robot_vector = state_all.block<3,1>(0,0);

            Eigen::MatrixXf line_robot(0,0);
            Eigen::MatrixXf R_line(0,0);
            Eigen::VectorXi line_idx(0);
            Eigen::MatrixXf endpoints_line_measure(0,0);

            line_match->performLineMatching(line_robot, R_line, line_idx, endpoints_line_measure, line_un_match, endpoints_un_match, R_un_match,\
        line_in_state, endpoints_in_state_use, line_test, endpoints_test, R_test, robot_vector);

            std::cout << "line_idx: \n" << line_idx << std::endl; 
            std::cout << "line_robot: \n" << line_robot << std::endl;


            ekf->measurement_update(state_all_new, Cov_k_k_new,\
        state_all,Cov,line_robot,R_line,line_idx);
            std::cout << "Measure_Cov:\n " << Cov_k_k_new(0,0) <<std::endl;
            Cov = Cov_k_k_new;
            state_all = state_all_new;
            state_all(2,0) = atan2(sin(state_all(2,0)),cos(state_all(2,0)));
            Eigen::MatrixXf endpoints_in_state_new;
            ekf->update_mached_line_endpoint(endpoints_in_state_new,\
        state_all, endpoints_in_state,line_idx,endpoints_line_measure);
            endpoints_in_state = endpoints_in_state_new;

        }

        //////// add the lines that are not in the state vector
        if (state_all.size() == 3 || line_un_match.size() > 0) {

            if (state_all.size() == 3) {
                line_add = line_para;
                R_add = R_mat;
                endpoints_add = endpoints_in_line;
            } else if (line_un_match.size() > 0) {
                line_add = line_un_match;
                R_add = R_un_match;
                endpoints_add = endpoints_un_match;
            }

            for (int i = 0; i < line_add.cols(); i++) {
                R_add_single = R_add.block(2*i,0,2,2);
                x_part << endpoints_add(0,i), endpoints_add(1,i);
                y_part << endpoints_add(2,i), endpoints_add(3,i);

                ekf->add_line_2_state(state_all_new,Cov_k_k_new,end_point_gl,\
          state_all,Cov,x_part,y_part,R_add_single);

                AddMatrixToCol(endpoints_in_state,end_point_gl);
                state_all = state_all_new;
                Cov = Cov_k_k_new;
            }

        }
        state_all(2,0) = atan2(sin(state_all(2,0)),cos(state_all(2,0)));
        //////// print endpoints
        for (int j = 0; j < endpoints_in_state.cols(); j++) {
            endFile<< endpoints_in_state(0,j)<<" "<<endpoints_in_state(1,j)\
            <<" "<<endpoints_in_state(2,j)<<" "<<endpoints_in_state(3,j)<< std::endl;
           end_update_file << endpoints_in_state(0,j)<<" "<<endpoints_in_state(1,j)\
            <<" "<<endpoints_in_state(2,j)<<" "<<endpoints_in_state(3,j) << " ";
	}
	if (endpoints_in_state.cols() > 0) {
	  end_update_file << std::endl;
	}

	std::cout << "dt: " << dt << std::endl;
	std::cout << "pose: x: " << state_all(0) << " y: " << state_all(1) << " theta: " << state_all(2) << std::endl << std::endl;
    
    poseFile << state_all(0) <<" "<< state_all(1) << std::endl;
    covFile << Cov(0,0)<<" "<<Cov(0,1)<<" "<<Cov(1,0)<<" "<<Cov(1,1)<< std::endl;

    loopTime += dt;

    if (loopTime > 1.0){      
      sick.lockDevice();
        std::vector<ArSensorReading> *r = sick.getRawReadingsAsVector();
        std::vector<ArSensorReading> readings(*r);
      sick.unlockDevice();
            
      for (int i=0; i<readings.size(); i++){
        if (readings[i].getRange() > 7000) continue;
        
        float fx = readings[i].getLocalX()/1000.0;
        float fy = readings[i].getLocalY()/1000.0;
        float newX = fx*cos(state_all(2)) - fy*sin(state_all(2));
        float newY = fx*sin(state_all(2)) + fy*cos(state_all(2));

        //Save scan to file
        laser_file << newX + state_all(0)/1000.0 << " " << newY + state_all(1)/1000.0 << std::endl;
      }
      loopTime = 0.0;
    }



    }

    robot.waitForRunExit();

    end_update_file.close();
    endFile.close();
    poseFile.close();
    covFile.close();
    laser_file.close();

    Aria::exit(0);
    return 0;

}

