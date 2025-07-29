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

    double proportional_to_rotation_loss(Mat3 &R);

    Mat34 add_noise_to_P(const Mat34 &P, double noise_multiplicator);

    class FocalLengths_DistributionAndLosses {
        
        // Each focal length has a dictionnary containing : Count, QR_Kcost, QR_Focalcost, QR_PPcost, QR_Skewcost and the IAC equivalents
        using FlLog = std::map<std::string, double>;
        using LossesMap = std::map<size_t, FlLog>;
    public:
        FocalLengths_DistributionAndLosses(double center, double stddev, bool use_LogNormal = false);
        void draw_random_fl(size_t num_cams); // Choose the focal lengths according to the distribution for the current reconstruction
        void add_loss(size_t idx, std::string which_loss, double loss); // Add correct loss corresponding to the idx-th image of the reconstrution
        void write_losses(std::string path); // Average losses and write them in a file for each focal length
        void write_losses_for_a_given_translation_range(std::ofstream &outFile, double translation_range, double Rank3Rate);
        void clear();

        std::vector<size_t> current_focal_lengths;

    private:
        double center;
        double stddev;
        bool use_LogNormal;

        LossesMap losses_map;
        std::normal_distribution<double> dist;
        std::lognormal_distribution<double> log_dist;
        std::mt19937 gen;
    };
}