function value = factorial_stirling ( n )

%% FACTORIAL_STIRLING computes Stirling's approximation to N!.
%
%  Definition:
%
%    N! = Product ( 1 <= I <= N ) I
%
%    Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )**N * E**(1/(12*N) )
%
%  Discussion:
%
%    This routine returns the raw approximation for all nonnegative
%    values of N.  If N is less than 0, the value is returned as 0,
%    and if N is 0, the value of 1 is returned.  In all other cases,
%    Stirling's formula is used.
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
%    Input, integer N, the argument of the function.
%
%    Output, real FACTORIAL_STIRLING, an approximation to N!.
%
  e_natural = 2.718281828459045;

  if ( n < 0 )

    value = 0.0;

  elseif ( n == 0 )

    value = 1.0;

  else

    value = sqrt ( 2.0 * pi * n ) * ( n / e_natural )^n * exp ( 1.0 / ( 12 * n ) );

  end
