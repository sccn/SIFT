function cdf = anglit_cdf ( x )

%% ANGLIT_CDF evaluates the Anglit CDF.
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
%    Output, real CDF, the value of the CDF.
%
  if ( x <  -0.25 * pi )
    cdf = 0.0;
  elseif ( x < 0.25 * pi )
    cdf = 0.5 - 0.5 * cos ( 2.0 * x + pi / 2.0 );
  else
    cdf = 1.0;
  end
