/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_observer_api.h
 *
 * Code generation for function '_coder_observer_api'
 *
 */

#ifndef _CODER_OBSERVER_API_H
#define _CODER_OBSERVER_API_H

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
extern void observer(real_T t, real_T theta, real_T varphi, real_T u, real_T q[2],
                     real_T dq[2]);
extern void observer_api(const mxArray * const prhs[4], int32_T nlhs, const
  mxArray *plhs[2]);
extern void observer_atexit(void);
extern void observer_initialize(void);
extern void observer_terminate(void);
extern void observer_xil_shutdown(void);
extern void observer_xil_terminate(void);
extern real_T test_controller(real_T theta, real_T varphi, real_T dtheta, real_T
  dvarphi);
extern void test_controller_api(const mxArray * const prhs[4], int32_T nlhs,
  const mxArray *plhs[1]);

#endif

/* End of code generation (_coder_observer_api.h) */
