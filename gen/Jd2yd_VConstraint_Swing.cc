/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:26 GMT-04:00
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
  double t1106;
  double t1107;
  double t1108;
  double t1120;
  double t1121;
  double t1125;
  double t1126;
  double t1127;
  double t1128;
  double t1129;
  double t1133;
  double t1134;
  double t1138;
  double t1139;
  double t1122;
  double t1123;
  double t1124;
  double t1130;
  double t1131;
  double t1132;
  double t1135;
  double t1136;
  double t1137;
  double t1140;
  double t1141;
  double t1142;
  double t1143;
  double t1147;
  double t1148;
  double t1149;
  double t1150;
  double t1151;
  double t1152;
  double t1153;
  double t1154;
  double t1155;
  double t1156;
  double t1157;
  double t1158;
  double t1159;
  double t1144;
  double t1145;
  double t1146;
  double t1160;
  double t1161;
  double t1162;
  double t1163;
  double t1172;
  double t1164;
  double t1165;
  double t1166;
  double t1167;
  double t1168;
  double t1169;
  double t1170;
  double t1171;
  double t1173;
  double t1174;
  double t1175;
  double t1176;
  double t1177;
  double t1178;
  double t1179;
  double t1180;
  double t1181;
  double t1182;
  double t1183;
  double t1184;
  double t1185;
  double t1186;
  t1106 = -1.*var4[0];
  t1107 = t1106 + var4[1];
  t1108 = Power(t1107,-5);
  t1120 = t1106 + var1[2] + var1[3];
  t1121 = Power(t1120,3);
  t1125 = Power(t1107,-4);
  t1126 = Power(t1120,2);
  t1127 = 1/t1107;
  t1128 = -1.*t1127*t1120;
  t1129 = 1. + t1128;
  t1133 = Power(t1107,-3);
  t1134 = Power(t1129,2);
  t1138 = Power(t1107,-2);
  t1139 = Power(t1129,3);
  t1122 = 20.*var3[6]*t1108*t1121;
  t1123 = -40.*var3[8]*t1108*t1121;
  t1124 = 20.*var3[10]*t1108*t1121;
  t1130 = 60.*var3[4]*t1125*t1126*t1129;
  t1131 = -120.*var3[6]*t1125*t1126*t1129;
  t1132 = 60.*var3[8]*t1125*t1126*t1129;
  t1135 = 60.*var3[2]*t1133*t1120*t1134;
  t1136 = -120.*var3[4]*t1133*t1120*t1134;
  t1137 = 60.*var3[6]*t1133*t1120*t1134;
  t1140 = 20.*var3[0]*t1138*t1139;
  t1141 = -40.*var3[2]*t1138*t1139;
  t1142 = 20.*var3[4]*t1138*t1139;
  t1143 = t1122 + t1123 + t1124 + t1130 + t1131 + t1132 + t1135 + t1136 + t1137 + t1140 + t1141 + t1142;
  t1147 = 20.*var3[7]*t1108*t1121;
  t1148 = -40.*var3[9]*t1108*t1121;
  t1149 = 20.*var3[11]*t1108*t1121;
  t1150 = 60.*var3[5]*t1125*t1126*t1129;
  t1151 = -120.*var3[7]*t1125*t1126*t1129;
  t1152 = 60.*var3[9]*t1125*t1126*t1129;
  t1153 = 60.*var3[3]*t1133*t1120*t1134;
  t1154 = -120.*var3[5]*t1133*t1120*t1134;
  t1155 = 60.*var3[7]*t1133*t1120*t1134;
  t1156 = 20.*var3[1]*t1138*t1139;
  t1157 = -40.*var3[3]*t1138*t1139;
  t1158 = 20.*var3[5]*t1138*t1139;
  t1159 = t1147 + t1148 + t1149 + t1150 + t1151 + t1152 + t1153 + t1154 + t1155 + t1156 + t1157 + t1158;
  t1144 = var2[2]*t1143;
  t1145 = var2[3]*t1143;
  t1146 = t1144 + t1145;
  t1160 = var2[2]*t1159;
  t1161 = var2[3]*t1159;
  t1162 = t1160 + t1161;
  t1163 = Power(t1120,4);
  t1172 = Power(t1129,4);
  t1164 = -5.*var3[8]*t1108*t1163;
  t1165 = 5.*var3[10]*t1108*t1163;
  t1166 = -20.*var3[6]*t1125*t1121*t1129;
  t1167 = 20.*var3[8]*t1125*t1121*t1129;
  t1168 = -30.*var3[4]*t1133*t1126*t1134;
  t1169 = 30.*var3[6]*t1133*t1126*t1134;
  t1170 = -20.*var3[2]*t1138*t1120*t1139;
  t1171 = 20.*var3[4]*t1138*t1120*t1139;
  t1173 = -5.*var3[0]*t1127*t1172;
  t1174 = 5.*var3[2]*t1127*t1172;
  t1175 = t1164 + t1165 + t1166 + t1167 + t1168 + t1169 + t1170 + t1171 + t1173 + t1174;
  t1176 = -5.*var3[9]*t1108*t1163;
  t1177 = 5.*var3[11]*t1108*t1163;
  t1178 = -20.*var3[7]*t1125*t1121*t1129;
  t1179 = 20.*var3[9]*t1125*t1121*t1129;
  t1180 = -30.*var3[5]*t1133*t1126*t1134;
  t1181 = 30.*var3[7]*t1133*t1126*t1134;
  t1182 = -20.*var3[3]*t1138*t1120*t1139;
  t1183 = 20.*var3[5]*t1138*t1120*t1139;
  t1184 = -5.*var3[1]*t1127*t1172;
  t1185 = 5.*var3[3]*t1127*t1172;
  t1186 = t1176 + t1177 + t1178 + t1179 + t1180 + t1181 + t1182 + t1183 + t1184 + t1185;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=t1146;
  p_output1[5]=t1162;
  p_output1[6]=t1146;
  p_output1[7]=t1162;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t1175;
  p_output1[15]=t1186;
  p_output1[16]=t1175;
  p_output1[17]=t1186;
  p_output1[18]=0;
  p_output1[19]=0;
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
    ( !(mrows == 12 && ncols == 1) && 
      !(mrows == 1 && ncols == 12))) 
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 10, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3,var4);


}

#else // MATLAB_MEX_FILE

#include "Jd2yd_VConstraint_Swing.hh"

namespace SymFunction
{

void Jd2yd_VConstraint_Swing_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4);

}

}

#endif // MATLAB_MEX_FILE
