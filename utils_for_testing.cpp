#include "utils_for_testing.hpp"
#include <cmath>
#include <iostream>

namespace rootba_povar {
    Mat3 RotationAroundX(double angle) {
        double c, s;
        sincos(angle, &s, &c);
        Mat3 R;
        R << 1,  0,  0,
            0,  c, -s,
            0,  s,  c;
        return R;
    }

    Mat3 RotationAroundY(double angle) {
        double c, s;
        sincos(angle, &s, &c);
        Mat3 R;
        R <<  c, 0, s,
            0, 1, 0,
            -s, 0, c;
        return R;
    }

    Mat3 RotationAroundZ(double angle) {
        double c, s;
        sincos(angle, &s, &c);
        Mat3 R;
        R << c, -s,  0,
            s,  c,  0,
            0,  0,  1;
        return R;
    }
    
    Mat3 random_rotation(double angle_range) {
        Mat3 R = RotationAroundX(double(rand()) / RAND_MAX * angle_range)
            * RotationAroundY(double(rand()) / RAND_MAX * angle_range)
            * RotationAroundZ(double(rand()) / RAND_MAX * angle_range);
        return R;
    }

    Vec3 random_translation(double translation_range) {
        Vec3 t(double(rand()) / RAND_MAX * translation_range,
            double(rand()) / RAND_MAX * translation_range,
            double(rand()) / RAND_MAX* translation_range);
        return t;
    }
    Mat34 P_out_of_random_Rt(Mat3 K, double translation_range ){
        Mat3 R = random_rotation();
        Vec3 t = random_translation(translation_range);
        Mat34 P_metric;
        P_From_KRt(K, R, t, &P_metric);
        return P_metric;
    }
    void P_From_KRt(const Mat3 &K, const Mat3 &R, const Vec3 &t, Mat34 *P) {
        P->block<3, 3>(0, 0) = R;
        P->col(3) = t;
        (*P) = K * (*P);
    }

    double Mat3_distance(Mat3 &A, Mat3 &B){
        // Normalized w.r.t A so that we get a distance invariant to the scale of the matrices compared.
        double tot = 0;
        for (int i=0; i<3; ++i){
            for (int j=0; j<3; ++j){
                double diff = A(i,j) - B(i,j);
                tot += diff * diff;
            }
        }
        return sqrt(tot);
    }

    double Mat3_distance_normalized(Mat3 &A, Mat3 &B){
        double tot_diff=0;
        double tot_A=0;
        double diff;
        for (int i=0; i<3; ++i){
            for (int j=0; j<3; ++j){
                diff = A(i,j) - B(i,j);
                tot_diff += diff * diff;
                tot_A += A(i,j) * A(i,j);
            }
        }
        tot_diff /= tot_A;
        tot_diff = sqrt(tot_diff);
        return tot_diff;
    }
    
    void K_From_ImageOfTheAbsoluteConic(Mat4 &Q, Mat34 &P, Mat3 *K){
        Mat3 KKt = P * Q * P.transpose();
        KKt /= KKt(2,2);
        K_From_ImageOfTheDualAbsoluteQuadric(KKt, K);
    }
    void random_Mat4(Mat4& P, double scale){
        P << double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale,
        double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale,
        double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale,
        double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale, double(rand()) / RAND_MAX * scale;
    }

    FocalLengths_DistributionAndLosses::FocalLengths_DistributionAndLosses(double center, double stddev)
        : center(center), stddev(stddev) 
    {
        std::random_device rd;
        gen = std::mt19937(rd());
        dist = std::normal_distribution<double>(center, stddev);
    }

    void FocalLengths_DistributionAndLosses::add_loss(size_t idx, std::string which_loss, double loss){
            losses_map[current_focal_lengths[idx]][which_loss] += loss;
    }

    void FocalLengths_DistributionAndLosses::draw_random_fl(size_t num_cams){
        current_focal_lengths.resize(num_cams);
        for (int i = 0; i < num_cams; ++i) {
            double number;
            do {number = dist(gen);} while (number < 1);
            size_t val = static_cast<size_t>(std::round(number)); // Choose integer focal length
            current_focal_lengths[i] = val;
            losses_map[val]["Count"] += 1;
        }
    }

    void FocalLengths_DistributionAndLosses::write_losses(std::string path){
        std::ofstream outFile(path);

        if (!outFile){
            std::cerr << "Erros opening file for writing!" << std::endl;
        }

        outFile << "FocalLength\tQR_Kcost\tQR_Focalcost\tQR_PPcost\tQR_Skewcost\tIAC_Kcost\tIAC_Focalcost\tIAC_PPcost\tIAC_Skewcost\n";
        for (auto& [focal_length, losses_list] : losses_map) {
            outFile << focal_length << "\t"
                    << losses_list["QR_Kcost"] / losses_list["Count"] << "\t"
                    << losses_list["QR_Focalcost"] / losses_list["Count"] << "\t"
                    << losses_list["QR_PPcost"] / losses_list["Count"] << "\t"
                    << losses_list["QR_Skewcost"] / losses_list["Count"] <<"\t"
                    << losses_list["IAC_Kcost"] / losses_list["Count"] << "\t"
                    << losses_list["IAC_Focalcost"] / losses_list["Count"] << "\t"
                    << losses_list["IAC_PPcost"] / losses_list["Count"] << "\t"
                    << losses_list["IAC_Skewcost"] / losses_list["Count"] <<"\n";
        }
        outFile.close();
    }
}