/*
 * Automatically Generated from Mathematica.
 * Tue 1 Aug 2017 13:02:21 GMT-04:00
 */

#ifndef KINEMATICSP_STANCE_JOINT_HH
#define KINEMATICSP_STANCE_JOINT_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymExpression
{

  void kinematicsp_stance_joint_raw(double *p_output1, const double *var1);

  inline void kinematicsp_stance_joint(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    kinematicsp_stance_joint_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // KINEMATICSP_STANCE_JOINT_HH
