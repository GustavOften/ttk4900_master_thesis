//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  observer_initialize.cpp
//
//  Code generation for function 'observer_initialize'
//


// Include files
#include "observer_initialize.h"
#include "observer.h"
#include "observer_data.h"
#include "rt_nonfinite.h"
#include "test_controller.h"

// Function Definitions
void observer_initialize()
{
  rt_InitInfAndNaN();
  spline_phi_not_empty_init();
  observer_init();
  test_controller_init();
  isInitialized_observer = true;
}

// End of code generation (observer_initialize.cpp)
