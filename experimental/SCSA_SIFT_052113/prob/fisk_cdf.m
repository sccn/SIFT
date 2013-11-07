function cdf = fisk_cdf ( x, a, b, c )

%% FISK_CDF evaluates the Fisk CDF.
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
%    0.0D+00 < B,
%    0.0D+00 < C.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= a )
    cdf = 0.0;
  else
    cdf = 1.0 / ( 1.0 + ( b / ( x - a ) )^c );
  end
