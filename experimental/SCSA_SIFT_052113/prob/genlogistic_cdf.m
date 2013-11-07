function cdf = genlogistic_cdf ( x, a, b, c )

%% GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
%
%  Modified:
%
%    11 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the CDF.
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real CDF, the value of the CDF.
%
  y = ( x - a ) / b;

  cdf = 1.0 / ( 1.0 + exp ( - y ) )^c;
