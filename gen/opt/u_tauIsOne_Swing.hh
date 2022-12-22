/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:39:55 GMT-04:00
 */

#ifndef U_TAUISONE_SWING_HH
#define U_TAUISONE_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace Pattern[ThreeLink, Blank[system]]
{

  void u_tauIsOne_Swing_raw(double *p_output1, const double *var1);

  inline void u_tauIsOne_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    u_tauIsOne_Swing_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // U_TAUISONE_SWING_HH
