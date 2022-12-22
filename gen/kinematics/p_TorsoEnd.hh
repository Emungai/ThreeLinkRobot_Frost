/*
 * Automatically Generated from Mathematica.
 * Fri 10 Jun 2022 15:48:55 GMT-04:00
 */

#ifndef P_TORSOEND_HH
#define P_TORSOEND_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymExpression
{

  void p_TorsoEnd_raw(double *p_output1, const double *var1);

  inline void p_TorsoEnd(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    p_TorsoEnd_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // P_TORSOEND_HH
