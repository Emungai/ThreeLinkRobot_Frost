/*
 * Automatically Generated from Mathematica.
 * Fri 10 Jun 2022 15:48:51 GMT-04:00
 */

#ifndef P_SWING_JOINT_HH
#define P_SWING_JOINT_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymExpression
{

  void p_swing_joint_raw(double *p_output1, const double *var1);

  inline void p_swing_joint(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    p_swing_joint_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // P_SWING_JOINT_HH
