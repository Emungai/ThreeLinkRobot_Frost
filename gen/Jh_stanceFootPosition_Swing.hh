/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:23 GMT-04:00
 */

#ifndef JH_STANCEFOOTPOSITION_SWING_HH
#define JH_STANCEFOOTPOSITION_SWING_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Jh_stanceFootPosition_Swing_raw(double *p_output1, const double *var1);

  inline void Jh_stanceFootPosition_Swing(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 5, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 5);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Jh_stanceFootPosition_Swing_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // JH_STANCEFOOTPOSITION_SWING_HH
