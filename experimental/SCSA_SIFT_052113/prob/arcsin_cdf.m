function cdf = arcsin_cdf ( x, a )

%% ARCSIN_CDF evaluates the Arcsin CDF.
%
%  Modified:
%
%    01 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the CDF.
%
%    Input, real A, the parameter of the CDF.
%    A must be positive.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= -a )
    cdf = 0.0;
  elseif ( x < a )
    cdf = 0.5 + asin ( x / a ) / pi;
  else
    cdf = 1.0;
  end
