#include "libmv.hpp"
#include "testing_functions.hpp"
#include <iostream>


int main(){
//   rootba_povar::test_K_From_AbsoluteConic();
//   rootba_povar::test_with_metric_input();
  // rootba_povar::recoverK_IACvsQR(10,100,1);
  // rootba_povar::recoverK_IACvsQR(50,100,1);
  // rootba_povar::recoverK_IACvsQR(10,100,10);
  // rootba_povar::recoverK_IACvsQR(10,100,10);
  // rootba_povar::recoverK_IACvsQR(50,100,10);
  // rootba_povar::recoverK_IACvsQR(10,100,100);
  // rootba_povar::recoverK_IACvsQR(50,100,100);
  // rootba_povar::recoverK_IACvsQR(10,100,1000);
  // rootba_povar::recoverK_IACvsQR(50,100,1000);
  // rootba_povar::recoverK_IACvsQR(10,100,10000);
  // rootba_povar::recoverK_IACvsQR(50,100,10000);

rootba_povar::recoverK_detailed_losses(50, 50, 1);
rootba_povar::recoverK_detailed_losses(50, 50, 100);

  // rootba_povar::recoverK_varying_focal_length(10, 5000, 1);
  // rootba_povar::recoverK_varying_focal_length(50, 5000, 1);

  // rootba_povar::recoverK_varying_focal_length(10, 5000, 100);

  // rootba_povar::recoverK_varying_focal_length(50, 5000, 100);

  // rootba_povar::recoverK_varying_focal_length(10, 5000, 1000);
  // rootba_povar::recoverK_varying_focal_length(50, 5000, 1000);

  
  return 0;
}
