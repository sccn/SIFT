function cnk_log = binomial_coef_log ( n, k )

%% BINOMIAL_COEF_LOG computes the logarithm of the Binomial coefficient.
%
%  Formula:
%
%    CNK_LOG = LOG ( C(N,K) ) = LOG ( N! / ( K! * (N-K)! ) ).
%
%  Modified:
%
%    03 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, K, are the values of N and K.
%
%    Output, real CNK_LOG, the logarithm of C(N,K).
%
  cnk_log = factorial_log ( n ) - factorial_log ( k ) - factorial_log ( n - k );
