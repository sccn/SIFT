function cdf = normal_cdf ( x, a, b )

%% NORMAL_CDF evaluates the Normal CDF.
%
%  Discussion:
%
%    The Normal CDF is related to the Error Function ERF(X) by:
%
%      ERF ( X ) = 2 * NORMAL_CDF ( SQRT ( 2 ) * X ) - 1.0.
%
%  Modified:
%
%    17 September 2004
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

  cdf = normal_01_cdf ( y );
