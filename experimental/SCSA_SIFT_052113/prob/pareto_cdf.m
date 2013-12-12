function cdf = pareto_cdf ( x, a, b )

%% PARETO_CDF evaluates the Pareto CDF.
%
%  Modified:
%
%    19 September 2004
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
%    0.0 < A,
%    0.0 < B.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x < a )
    cdf = 0.0;
  else
    cdf = 1.0 - ( a / x )^b;
  end
