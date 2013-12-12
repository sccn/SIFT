function cdf = exponential_cdf ( x, a, b )

%% EXPONENTIAL_CDF evaluates the Exponential CDF.
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
%    Input, real X, the argument of the PDF.
%
%    Input, real A, B, the parameter of the PDF.
%    0.0 < B.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= a )
    cdf = 0.0;
  else
    cdf = 1.0 - exp ( ( a - x ) / b );
  end

