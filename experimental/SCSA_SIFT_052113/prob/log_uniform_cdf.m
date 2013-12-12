function cdf = log_uniform_cdf ( x, a, b )

%% LOG_UNIFORM_CDF evaluates the Log Uniform CDF.
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
%    Input, real X, the argument of the CDF.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= a )
    cdf = 0.0;
  elseif ( x < b )
    cdf = ( log ( x ) - log ( a ) ) / ( log ( b ) - log ( a ) );
  else
    cdf = 1.0;
  end
