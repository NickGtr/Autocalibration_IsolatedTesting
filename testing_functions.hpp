// Nicolas : Used to test the performance of autocalibration. Not templated : double is used

#pragma once

#include "libmv.hpp"

namespace rootba_povar {
    void test_K_From_AbsoluteConic();

    void test_with_metric_input();

    void recoverK_IACvsQR(int num_cams = 10, int num_reconstructions = 1, double translation_range = 1);

    void recoverK_detailed_losses(int num_cams = 10, int num_reconstructions = 10, double translation_range = 1);

    void recoverK_varying_focal_length(int num_cams = 10, int num_reconstructions = 100, double translation_range = 1);

    void recoverK_FIXED_IACvsQR_effect_of_translation_range(int num_cams, int num_reconstructions, int equivalent_focal_length);

    void recoverK_VARYING_IACvsQR_effect_of_translation_range(int num_cams, int num_reconstructions);
}