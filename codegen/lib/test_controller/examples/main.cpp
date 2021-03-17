//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  main.cpp
//
//  Code generation for function 'main'
//


//***********************************************************************
// This automatically generated example C++ main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************

// Include files
#include "main.h"
#include "rt_nonfinite.h"
#include "test_controller.h"
#include "test_controller_terminate.h"

// Function Declarations
static double argInit_real_T();
static void main_test_controller();

// Function Definitions
static double argInit_real_T()
{
  return 0.0;
}

static void main_test_controller()
{
  double theta_tmp;

  // Initialize function 'test_controller' input arguments.
  theta_tmp = argInit_real_T();

  // Call the entry-point 'test_controller'.
  theta_tmp = test_controller(theta_tmp, theta_tmp, theta_tmp, theta_tmp);
}

int main(int, const char * const [])
{
  // The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. 
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_test_controller();

  // Terminate the application.
  // You do not need to do this more than one time.
  test_controller_terminate();
  return 0;
}

// End of code generation (main.cpp)
