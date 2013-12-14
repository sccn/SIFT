function value = gamma_log_int ( n )

%% GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
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
%    Input, integer N, the argument of the logarithm of the Gamma function.
%    0 < N.
%
%    Output, real VALUE, the logarithm of
%    the Gamma function of N.
%
  if ( n <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GAMMA_LOG_INT - Fatal error!\n' );
    fprintf ( 1, '  Illegal input value of N = %d\n', n );
    fprintf ( 1, '  But N must be strictly positive.\n' );
    error ( 'GAMMA_LOG_INT - Fatal error!' );
  end

  value = gamma_log ( n );
