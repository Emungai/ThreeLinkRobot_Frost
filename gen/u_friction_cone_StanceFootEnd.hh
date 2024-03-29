/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 15:57:03 GMT-04:00
 */

#ifndef U_FRICTION_CONE_STANCEFOOTEND_HH
#define U_FRICTION_CONE_STANCEFOOTEND_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void u_friction_cone_StanceFootEnd_raw(double *p_output1, const double *var1,const double *var2);

  inline void u_friction_cone_StanceFootEnd(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 3, 1);
    assert_size_matrix(var2, 1, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 5, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    u_friction_cone_StanceFootEnd_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // U_FRICTION_CONE_STANCEFOOTEND_HH
