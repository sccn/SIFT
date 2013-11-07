function cdf = logistic_cdf ( x, a, b )

%% LOGISTIC_CDF evaluates the Logistic CDF.
%
%  Modified:
%
%    13 September 2004
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
  cdf = 1.0 / ( 1.0 + exp ( ( a - x ) / b ) );
