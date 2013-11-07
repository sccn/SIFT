function cdf = angle_cdf ( x, n )

%% ANGLE_CDF evaluates the Angle CDF.
%
%  Modified:
%
%    09 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Reuven Rubinstein,
%    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
%    Wiley, 1986.
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%
%    Input, integer N, the spatial dimension.
%    N must be at least 2.
%
%    Output, real CDF, the value of the CDF.
%
  if ( n < 2 )
    fprintf ( 1, '\n' );
    fprintf ( 1, '\ANGLE_CDF - Fatal error!\n' );
    fprintf ( 1, '\  N must be at least 2.\n' );
    fprintf ( 1, '\  The input value of N = %d\n', n );
    error ( 'ANGLE_PDF - Fatal error!' );
  end

  if ( x <= 0.0 )
    cdf = 0.0;
  elseif ( pi <= x )
    cdf = 1.0;
  elseif ( n == 2 )
    cdf = x / pi;
  else
    cdf = sin_power_int ( 0.0, x, n - 2 ) * gamma ( n / 2.0 ) ...
      / ( sqrt ( pi ) * gamma ( ( n - 1 ) / 2.0 ) );
  end
