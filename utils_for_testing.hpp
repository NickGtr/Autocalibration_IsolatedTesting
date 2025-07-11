// Nicolas : tools for testing in testing_functions

#pragma once

#include "libmv.hpp"
#include <random>
#include <vector>
#include <map>
#include <string>
#include <fstream>

namespace rootba_povar {
    Mat3 RotationAroundX(double angle);
    Mat3 RotationAroundY(double angle);
    Mat3 RotationAroundZ(double angle);
    Mat3 random_rotation(double angle_range = 3);
    Vec3 random_translation(double translation_range = 1);
    Mat34 P_out_of_random_Rt(Mat3 K, double translation_range = 1);

    void P_From_KRt(const Mat3 &K, const Mat3 &R, const Vec3 &t, Mat34 *P);

    void K_From_ImageOfTheAbsoluteConic(Mat4 &Q, Mat34 &P, Mat3 *K);

    double Mat3_distance(Mat3 &A, Mat3 &B);

    void random_Mat4(Mat4& P, double scale = 10);

    class FocalLengths_DistributionAndLosses {
        
        // Each focal length has a dictionnary containing : Count, QR_Kcost, QR_Focalcost, QR_PPcost, QR_Skewcost and the IAC equivalents
        using FlLog = std::map<std::string, double>;
        using LossesMap = std::map<size_t, FlLog>;
    public:
        FocalLengths_DistributionAndLosses(double center, double stddev);
        void draw_random_fl(size_t num_cams); // Choose the focal lengths according to the distribution for the current reconstruction
        void add_loss(size_t idx, std::string which_loss, double loss); // Add correct loss corresponding to the idx-th image of the reconstrution
        void write_losses(std::string path); // Average losses and write them in a file

        std::vector<size_t> current_focal_lengths;

    private:
        double center;
        double stddev;

        LossesMap losses_map;
        std::normal_distribution<double> dist;
        std::mt19937 gen;
    };
}