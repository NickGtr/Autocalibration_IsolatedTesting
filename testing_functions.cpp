#include "testing_functions.hpp"
#include "utils_for_testing.hpp"

#include <iostream>
#include <cmath>
#include <sstream>

namespace rootba_povar {
    void test_K_From_AbsoluteConic(){
        Mat3 K, Kp;
            K << 10,  1, 30,
                0, 20, 40,
                0,  0,  1;

        Mat3 w = (K * K.transpose()).inverse();
        K_From_AbsoluteConic(w, &Kp);

        std::cout << Kp << std::endl;
        double cost = Mat3_distance(K, Kp);
        std::cout << cost << std::endl;
    }
    
    void test_with_metric_input(){
        double width = 1000, height = 800;
        Mat3 K;
        K << width,     0,  width / 2, // 1000x800 image with 35mm equiv focal length.
                0, width, height / 2,
                0,     0,          1;

        AutoCalibrationLinear<double> a;

        // Add cameras with random rotation and translation.
        for (int i = 0; i < 3; ++i) {
            a.AddProjection(P_out_of_random_Rt(K), width, height);
        }
  
        // Compute metric update transformation.
        Mat4 H = a.MetricTransformation();

        // Since the input was metric, the transformation should be a similarity.
        // The 3x3 submatrix should proportional to an orthonormal matrix.
        Mat3 R = H.block<3, 3>(0, 0);
        Mat3 RRt = R * R.transpose();
        std::cout << "The following matrix should be proportional to identity :" << std::endl;
        std::cout << RRt << std::endl;
    }

    void recoverK_IACvsQR(int num_cams, int num_reconstructions, double translation_range){
        // A FIXED focal_length for each reconstruction.
        // We try to recover K from a random projective distortions by comparing IAC and QR to get K from the metric upgrade.

        // By using the real dimensions of a full-frame sensor, the focal length is equal to the equivalent focal length
        const double width = 36, height = 24;

        std::stringstream path;
        path << "/Data/gautt/bundle_adjustment/autocalibration_project/results/QRvsIAC_"
            << num_cams << "Cams_"
            << num_reconstructions << "Reconstr_"
            << translation_range << "Transl.txt";
        // path << "/Data/gautt/bundle_adjustment/autocalibration_project/results/QRvsIAC_focal_length_" << num_cams << "cams.txt";
        std::cout << "writing : " << path.str() << std::endl;
        std::ofstream outFile(path.str());

        if (!outFile){
            std::cerr << "Erros opening file for writing!" << std::endl;
        }

        outFile << "FocalLength\tQRloss\tIACloss\n";

        for (size_t equivalent_focal_length = 10; equivalent_focal_length < 300; ++equivalent_focal_length){
            Mat3 K;
            K << equivalent_focal_length, 0, width / 2,
                0, equivalent_focal_length, height / 2,
                0, 0, 1;

            double average_QR_cost = 0, average_iac_cost = 0;

            for (size_t reconstruction_counter = 0; reconstruction_counter < num_reconstructions; ++reconstruction_counter){
                Mat4 H_real;
                random_Mat4(H_real);
                AutoCalibrationLinear<double> a;

                // Add cameras with random rotation and translation.
                Mat34 Ps[num_cams];
                for (int i = 0; i < num_cams; ++i) {
                    Ps[i] = P_out_of_random_Rt(K, translation_range) * H_real.inverse();  // Distort cameras.
                    a.AddProjection(Ps[i], width, height);
                }
    
                // Compute metric update transformation. Also recover the absolute quadric to try to recover K from it via IAC and compare with QR 
                Mat4 Q;
                Mat4 H_computed = a.MetricTransformation(&Q);
        
                for (int i = 0; i < num_cams; ++i) {
                    Mat34 P_metric = Ps[i] * H_computed;  // Undistort cameras.
                    Mat3 K_from_QR, K_from_iac, R;
                    Vec3 t;

                    KRt_From_P(P_metric, &K_from_QR, &R, &t);
                    K_From_ImageOfTheAbsoluteConic(Q, Ps[i], &K_from_iac);

                    average_QR_cost += Mat3_distance(K, K_from_QR); 
                    average_iac_cost += Mat3_distance(K, K_from_iac);
                }
            }
            average_QR_cost /= (num_reconstructions * num_cams);
            average_iac_cost /= (num_reconstructions * num_cams);

            outFile << equivalent_focal_length << "\t" << average_QR_cost << "\t" << average_iac_cost << "\n";
        }
        outFile.close();
    }

    void recoverK_detailed_losses(int num_cams, int num_reconstructions, double translation_range){
        // A FIXED focal length for each reconstruction.
        // Outputs detailed losses (in global, focal_length, principal point and skew) when recovering K through QR and IAC

        // By using the real dimensions of a full-frame sensor, the focal length is equal to the equivalent focal length
        const double width = 36, height = 24;
        // Max range of camera_translation when picking a translation at random

        std::stringstream path;
        path << "/Data/gautt/bundle_adjustment/autocalibration_project/results/RecoverK&R_detailed_"
             << num_cams << "Cams_"
             << num_reconstructions << "Reconstr_"
             << translation_range << "Transl.txt";

        std::ofstream outFile(path.str());

        if (!outFile){
            std::cerr << "Errors opening file for writing!" << std::endl;
        }

        outFile << "FocalLength\tQR_Kcost\tQR_Focalcost\tQR_PPcost\tQR_Skewcost\tQR_RRtcost\tIAC_Kcost\tIAC_Focalcost\tIAC_PPcost\tIAC_Skewcost\tIAC_RRtcost\tRank3Rate\n";

        for (size_t equivalent_focal_length = 10; equivalent_focal_length < 300; ++equivalent_focal_length){
            Mat3 K;
            K << equivalent_focal_length, 0, width / 2,
                0, equivalent_focal_length, height / 2,
                0, 0, 1;

            double  averageQR_K_cost = 0, averageQR_focal_cost = 0, averageQR_pp_cost = 0 ,averageQR_skew_cost = 0, averageQR_RRt_cost = 0,
                    averageIAC_K_cost = 0, averageIAC_focal_cost = 0, averageIAC_pp_cost = 0 ,averageIAC_skew_cost = 0, averageIAC_RRt_cost = 0,
                    rank3_rate = 0;
            
            for (size_t reconstruction_counter = 0; reconstruction_counter < num_reconstructions; ++reconstruction_counter){
                Mat4 H_real;
                random_Mat4(H_real);
                AutoCalibrationLinear<double> a;

                // Add cameras with random rotation and translation.
                Mat34 Ps[num_cams];
                for (int i = 0; i < num_cams; ++i) {                   
                    Ps[i] = P_out_of_random_Rt(K, translation_range) * H_real.inverse();  // Distort cameras.
                    a.AddProjection(Ps[i], width, height);
                }
    
                // Compute metric update transformation. Also recover the absolute quadric to try to recover K from it via IAC and compare with QR 
                Mat4 Q;
                int rank = 0;
                Mat4 H_computed = a.MetricTransformation(&Q, &rank);
                rank3_rate += (rank == 3);
        
                for (int i = 0; i < num_cams; ++i) {
                    Mat34 P_metric = Ps[i] * H_computed;  // Undistort cameras.
                    Mat3 K_from_QR, K_from_iac, R_from_QR, R_from_iac;
                    Vec3 t_from_QR;

                    KRt_From_P(P_metric, &K_from_QR, &R_from_QR, &t_from_QR);
                    K_From_ImageOfTheAbsoluteConic(Q, Ps[i], &K_from_iac);

                    R_from_iac = K_from_iac.inverse() * P_metric.block<3, 3>(0, 0);

                    averageQR_K_cost += Mat3_distance(K, K_from_QR);
                    averageQR_focal_cost += sqrt((K(1,1) - K_from_QR(1,1)) * (K(1,1) - K_from_QR(1,1)) + (K(0,0) - K_from_QR(0,0)) * (K(0,0) - K_from_QR(0,0)));
                    averageQR_pp_cost += sqrt((K(0,2) - K_from_QR(0,2)) * (K(0,2) - K_from_QR(0,2)) + (K(1,2) - K_from_QR(1,2)) * (K(1,2) - K_from_QR(1,2)));
                    averageQR_skew_cost += (K(0,1) - K_from_QR(0,1)) * (K(0,1) - K_from_QR(0,1));
                    averageQR_RRt_cost += proportional_to_rotation_loss(R_from_QR);

                    averageIAC_K_cost += Mat3_distance(K, K_from_iac);
                    averageIAC_focal_cost += sqrt((K(1,1) - K_from_iac(1,1)) * (K(1,1) - K_from_iac(1,1)) + (K(0,0) - K_from_iac(0,0)) * (K(0,0) - K_from_iac(0,0)));
                    averageIAC_pp_cost += sqrt((K(0,2) - K_from_iac(0,2)) * (K(0,2) - K_from_iac(0,2)) + (K(1,2) - K_from_iac(1,2)) * (K(1,2) - K_from_iac(1,2)));
                    averageIAC_skew_cost += (K(0,1) - K_from_iac(0,1)) * (K(0,1) - K_from_iac(0,1));
                    averageIAC_RRt_cost += proportional_to_rotation_loss(R_from_iac);
                }
            }
            averageQR_K_cost /= (num_reconstructions * num_cams);
            averageQR_focal_cost /= (num_reconstructions * num_cams);
            averageQR_pp_cost /= (num_reconstructions * num_cams);
            averageQR_skew_cost /= (num_reconstructions * num_cams);
            averageQR_RRt_cost /= (num_reconstructions * num_cams);

            averageIAC_K_cost /= (num_reconstructions * num_cams);
            averageIAC_focal_cost /= (num_reconstructions * num_cams);
            averageIAC_pp_cost /= (num_reconstructions * num_cams);
            averageIAC_skew_cost /= (num_reconstructions * num_cams);
            averageIAC_RRt_cost /= (num_reconstructions * num_cams);
            rank3_rate /= num_reconstructions;

            outFile << equivalent_focal_length << "\t"
                    << averageQR_K_cost << "\t"
                    << averageQR_focal_cost << "\t"
                    << averageQR_pp_cost << "\t"
                    << averageQR_skew_cost <<"\t"
                    << averageQR_RRt_cost << "\t"
                    << averageIAC_K_cost << "\t"
                    << averageIAC_focal_cost << "\t"
                    << averageIAC_pp_cost << "\t"
                    << averageIAC_skew_cost << "\t"
                    << averageIAC_RRt_cost << "\t"
                    << rank3_rate << "\n";
        }
        outFile.close();
    }

    void recoverK_varying_focal_length(int num_cams, int num_reconstructions, double translation_range){
        // VARYING focal lengths for each reconstruction
        // Outputs detailed losses (in global, focal_length, principal point and skew) when recovering K through QR and IAC.

        // Here we pick the distribution of focal length of the cameras (Normal or log-normal))
        const double center = 3.7, stddev = 0.6;
        const bool use_log_normal = true; // If false, use normal distribution.
        // By using the real dimensions of a full-frame sensor, the focal length is equal to the equivalent focal length
        const double width = 36, height = 24; 
    
        Mat3 K;
        Mat4 H_real;
        FocalLengths_DistributionAndLosses fl(center, stddev, use_log_normal);

        for (size_t reconstruction_counter = 0; reconstruction_counter < num_reconstructions; ++reconstruction_counter){
            Mat4 H_real;
            random_Mat4(H_real);
            AutoCalibrationLinear<double> a;
            fl.draw_random_fl(num_cams); //Focal_length used in reconstruction vary.

            Mat34 Ps[num_cams];
            Mat3 Ks[num_cams];
            for (size_t i = 0; i < num_cams; ++i){
                Ks[i] << fl.current_focal_lengths[i], 0, width / 2,
                            0, fl.current_focal_lengths[i], height / 2,
                            0, 0, 1;

                Ps[i] = P_out_of_random_Rt(Ks[i], translation_range) * H_real.inverse(); // Distort cameras
                a.AddProjection(Ps[i], width, height);
            }

            Mat4 Q;
            Mat4 H_computed = a.MetricTransformation(&Q);

            for (size_t i = 0; i < num_cams; ++i) {
                Mat34 P_metric = Ps[i] * H_computed;  // Undistort cameras.
                Mat3 K_from_QR, K_from_iac, R;
                Vec3 t;

                KRt_From_P(P_metric, &K_from_QR, &R, &t);
                K_From_ImageOfTheAbsoluteConic(Q, Ps[i], &K_from_iac);

                fl.add_loss(i, "QR_Kcost", Mat3_distance(Ks[i], K_from_QR));
                fl.add_loss(i, "QR_Focalcost", sqrt((Ks[i](1,1) - K_from_QR(1,1)) * (Ks[i](1,1) - K_from_QR(1,1)) + (Ks[i](0,0) - K_from_QR(0,0)) * (Ks[i](0,0) - K_from_QR(0,0))));
                fl.add_loss(i, "QR_PPcost", sqrt((Ks[i](0,2) - K_from_QR(0,2)) * (Ks[i](0,2) - K_from_QR(0,2)) + (Ks[i](1,2) - K_from_QR(1,2)) * (Ks[i](1,2) - K_from_QR(1,2))));
                fl.add_loss(i, "QR_Skewcost", (Ks[i](0,1) - K_from_QR(0,1)) * (Ks[i](0,1) - K_from_QR(0,1)));

                fl.add_loss(i, "IAC_Kcost", Mat3_distance(Ks[i], K_from_iac));
                fl.add_loss(i, "IAC_Focalcost", sqrt((Ks[i](1,1) - K_from_iac(1,1)) * (Ks[i](1,1) - K_from_iac(1,1)) + (Ks[i](0,0) - K_from_iac(0,0)) * (Ks[i](0,0) - K_from_iac(0,0))));
                fl.add_loss(i, "IAC_PPcost", sqrt((Ks[i](0,2) - K_from_iac(0,2)) * (Ks[i](0,2) - K_from_iac(0,2)) + (Ks[i](1,2) - K_from_iac(1,2)) * (Ks[i](1,2) - K_from_iac(1,2))));
                fl.add_loss(i, "IAC_Skewcost", (Ks[i](0,1) - K_from_iac(0,1)) * (Ks[i](0,1) - K_from_iac(0,1)));
            }
        }
        std::stringstream path;
        path << "/Data/gautt/bundle_adjustment/autocalibration_project/results/cRecoverK_varying_detailed_"
             << num_cams << "Cams_"
             << num_reconstructions << "Reconstr_"
             << translation_range << "Transl"
             << "_CamDistLogN(" << center << "," << stddev << ").txt";

        fl.write_losses(path.str());
    }
}