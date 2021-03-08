/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_test_controller_api.h
 *
 * Code generation for function '_coder_test_controller_api'
 *
 */

#ifndef _CODER_TEST_CONTROLLER_API_H
#define _CODER_TEST_CONTROLLER_API_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T test_controller(real_T theta, real_T varphi, real_T dtheta, real_T
  dvarphi);
extern void test_controller_api(const mxArray * const prhs[4], int32_T nlhs,
  const mxArray *plhs[1]);
extern void test_controller_atexit(void);
extern void test_controller_initialize(void);
extern void test_controller_terminate(void);
extern void test_controller_xil_shutdown(void);
extern void test_controller_xil_terminate(void);

#endif

/* End of code generation (_coder_test_controller_api.h) */
