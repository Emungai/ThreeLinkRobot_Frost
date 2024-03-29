/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:27 GMT-04:00
 */

#ifndef AVCONSTRAINTCONT_SWING_HH
#define AVCONSTRAINTCONT_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void aVConstraintCont_Swing_raw(double *p_output1, const double *var1,const double *var2);

  inline void aVConstraintCont_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 1, 12);
    assert_size_matrix(var2, 1, 12);

	
    // - Outputs
    assert_size_matrix(p_output1, 12, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    aVConstraintCont_Swing_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // AVCONSTRAINTCONT_SWING_HH
