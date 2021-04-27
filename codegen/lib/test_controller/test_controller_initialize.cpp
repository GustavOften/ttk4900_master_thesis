//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  test_controller_initialize.cpp
//
//  Code generation for function 'test_controller_initialize'
//


// Include files
#include "test_controller_initialize.h"
#include "rt_nonfinite.h"
#include "test_controller.h"
#include "test_controller_data.h"

// Function Definitions
void test_controller_initialize()
{
  rt_InitInfAndNaN();
  spline_phi_not_empty_init();
  test_controller_init();
  isInitialized_test_controller = true;
}

// End of code generation (test_controller_initialize.cpp)
