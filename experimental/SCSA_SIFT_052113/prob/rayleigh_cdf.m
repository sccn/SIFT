function cdf = rayleigh_cdf ( x, a )

%% RAYLEIGH_CDF evaluates the Rayleigh CDF.
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
%    0.0 <= X.
%
%    Input, real A, the parameter of the PDF.
%    0.0 < A.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x < 0.0 )
    cdf = 0.0;
  else
    cdf = 1.0 - exp ( - x * x / ( 2.0 * a * a ) );
  end
