function cdf = sech_cdf ( x, a, b )

%% SECH_CDF evaluates the Hyperbolic Secant CDF.
%
%  Modified:
%
%    20 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%
%    Input, real A, B, the parameter of the PDF.
%    0.0 < B.
%
%    Output, real CDF, the value of the CDF.
%
  y = ( x - a ) / b;

  cdf = 2.0 * atan ( exp ( y ) ) / pi;
