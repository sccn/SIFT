function cdf = reciprocal_cdf ( x, a, b )

%% RECIPROCAL_CDF evaluates the Reciprocal CDF.
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
%    Input, real A, B, the parameters of the PDF.
%    0.0D+00 < A <= B.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= 0.0 )

    cdf = 0.0;

  elseif ( 0.0 < x )

    cdf = log ( a / x ) / log ( a / b );

  end
