function value = beta ( a, b )

%% BETA returns the value of the Beta function.
%
%  Formula:
%
%    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
%              = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
%
%  Modified:
%
%    02 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the function.
%    0.0D+00 < A,
%    0.0D+00 < B.
%
%    Output, real BETA, the value of the function.
%
  if ( a <= 0.0 | b <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'BETA - Fatal error!\n' );
    fprintf ( 1, '  Both A and B must be greater than 0.\n' );
    error ( 'BETA - Fatal error!' );
  end

  value = exp ( gamma_log ( a ) + gamma_log ( b ) - gamma_log ( a + b ) );
