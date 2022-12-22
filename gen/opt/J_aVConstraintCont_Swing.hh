/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:28 GMT-04:00
 */

#ifndef J_AVCONSTRAINTCONT_SWING_HH
#define J_AVCONSTRAINTCONT_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void J_aVConstraintCont_Swing_raw(double *p_output1, const double *var1,const double *var2);

  inline void J_aVConstraintCont_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 1, 12);
    assert_size_matrix(var2, 1, 12);

	
    // - Outputs
    assert_size_matrix(p_output1, 24, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    J_aVConstraintCont_Swing_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // J_AVCONSTRAINTCONT_SWING_HH
