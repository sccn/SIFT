function cdf = extreme_values_cdf ( x, a, b )

%% EXTREME_VALUES_CDF evaluates the Extreme Values CDF.
%
%  Modified:
%
%    09 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the CDF.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real CDF, the value of the CDF.
%
  y = ( x - a ) / b;

  cdf = exp ( - exp ( - y ) );
