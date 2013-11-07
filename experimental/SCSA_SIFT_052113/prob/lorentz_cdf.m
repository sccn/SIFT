function cdf = lorentz_cdf ( x )

%% LORENTZ_CDF evaluates the Lorentz CDF.
%
%  Modified:
%
%    01 February 1999
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the CDF.
%
%    Output, real CDF, the value of the CDF.
%
  cdf = 0.5 + atan ( x ) / pi;
