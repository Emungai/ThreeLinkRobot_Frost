/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:45:15 GMT-04:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;


#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3,const double *var4)
{
  double t8734;
  double t8737;
  double t8735;
  double t8738;
  double t8736;
  double t8739;
  double t8740;
  double t8744;
  double t8746;
  double t8745;
  double t8747;
  double t8748;
  double t8756;
  double t8757;
  double t8758;
  double t8764;
  double t8765;
  double t8766;
  double t8731;
  double t8733;
  double t8752;
  double t8753;
  double t8754;
  double t8742;
  double t8743;
  double t8760;
  double t8761;
  double t8762;
  double t8770;
  double t8771;
  double t8755;
  double t8759;
  double t8763;
  double t8767;
  double t8768;
  double t8750;
  double t8751;
  double t8772;
  double t8773;
  double t8774;
  double t8775;
  double t8778;
  double t8779;
  double t8782;
  double t8783;
  double t8786;
  double t8787;
  double t8789;
  double t8790;
  double t8791;
  double t8792;
  double t8793;
  double t8795;
  double t8796;
  double t8797;
  double t8798;
  double t8799;
  double t8814;
  double t8815;
  double t8816;
  double t8817;
  double t8788;
  double t8794;
  double t8800;
  double t8801;
  double t8820;
  double t8821;
  double t8822;
  double t8823;
  double t8824;
  t8734 = Cos(var2[2]);
  t8737 = Sin(var2[2]);
  t8735 = Cos(var2[3]);
  t8738 = Sin(var2[3]);
  t8736 = t8734*t8735;
  t8739 = -1.*t8737*t8738;
  t8740 = t8736 + t8739;
  t8744 = Cos(var2[4]);
  t8746 = Sin(var2[4]);
  t8745 = t8734*t8744;
  t8747 = -1.*t8737*t8746;
  t8748 = t8745 + t8747;
  t8756 = t8735*t8737;
  t8757 = t8734*t8738;
  t8758 = t8756 + t8757;
  t8764 = t8744*t8737;
  t8765 = t8734*t8746;
  t8766 = t8764 + t8765;
  t8731 = -1.*var1[4];
  t8733 = var3[3] + t8731;
  t8752 = -1.*t8735*t8737;
  t8753 = -1.*t8734*t8738;
  t8754 = t8752 + t8753;
  t8742 = -1.*var1[3];
  t8743 = var3[4] + t8742;
  t8760 = -1.*t8744*t8737;
  t8761 = -1.*t8734*t8746;
  t8762 = t8760 + t8761;
  t8770 = -1.*var1[0];
  t8771 = var3[0] + t8770;
  t8755 = 5.*t8754*t8740;
  t8759 = 5.*t8758*t8740;
  t8763 = 5.*t8762*t8748;
  t8767 = 5.*t8766*t8748;
  t8768 = t8755 + t8759 + t8763 + t8767;
  t8750 = -1.*var1[1];
  t8751 = var3[1] + t8750;
  t8772 = Power(t8734,2);
  t8773 = 25.*t8772;
  t8774 = Power(t8737,2);
  t8775 = 25.*t8774;
  t8778 = Power(t8740,2);
  t8779 = 5.*t8778;
  t8782 = Power(t8748,2);
  t8783 = 5.*t8782;
  t8786 = -1.*var1[2];
  t8787 = var3[2] + t8786;
  t8789 = Power(t8735,2);
  t8790 = -0.5*t8789;
  t8791 = Power(t8738,2);
  t8792 = -0.5*t8791;
  t8793 = t8790 + t8792;
  t8795 = Power(t8744,2);
  t8796 = -0.5*t8795;
  t8797 = Power(t8746,2);
  t8798 = -0.5*t8797;
  t8799 = t8796 + t8798;
  t8814 = -10.*t8737;
  t8815 = 5.*t8754*t8793;
  t8816 = 5.*t8762*t8799;
  t8817 = t8814 + t8815 + t8816;
  t8788 = 10.*t8734;
  t8794 = 5.*t8740*t8793;
  t8800 = 5.*t8748*t8799;
  t8801 = t8788 + t8794 + t8800;
  t8820 = -1.*var4[1]*t8758;
  t8821 = -1.*t8734*t8735;
  t8822 = t8737*t8738;
  t8823 = t8821 + t8822;
  t8824 = -1.*var4[0]*t8823;
  p_output1[0]=-2.5*t8733*t8740 - 2.5*t8743*t8748 + t8751*t8768 + t8771*(5.*Power(t8758,2) + 5.*Power(t8766,2) + t8773 + t8775 + t8779 + t8783) + t8787*t8801 - 1.*var4[0];
  p_output1[1]=-2.5*t8733*t8754 - 2.5*t8743*t8762 + t8768*t8771 + t8751*(5.*Power(t8754,2) + 5.*Power(t8762,2) + t8773 + t8775 + t8779 + t8783) + t8787*t8817 - 1.*var4[1];
  p_output1[2]=-2.5*t8733*t8793 - 2.5*t8743*t8799 + t8787*(4.47 + 5.*Power(t8793,2) + 5.*Power(t8799,2)) + t8771*t8801 + t8751*t8817 + t8820 + t8824;
  p_output1[3]=1.25*t8733 - 2.5*t8751*t8754 - 2.5*t8740*t8771 - 2.5*t8787*t8793 + t8820 + t8824;
  p_output1[4]=1.25*t8743 - 2.5*t8751*t8762 - 2.5*t8748*t8771 - 2.5*t8787*t8799;
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1,*var2,*var3,*var4;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 4)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Four input(s) required (var1,var2,var3,var4).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var3 is wrong.");
    }
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
    ( !(mrows == 2 && ncols == 1) && 
      !(mrows == 1 && ncols == 2))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var4 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
  var4 = mxGetPr(prhs[3]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 5, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3,var4);


}

#else // MATLAB_MEX_FILE

#include "dxDiscreteMapSwingLegImpact.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void dxDiscreteMapSwingLegImpact_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4);

}

}

#endif // MATLAB_MEX_FILE
