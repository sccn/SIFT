function value = euler_constant ( dummy )

%% EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
%
%  Discussion:
%
%    The Euler-Mascheroni constant is often denoted by a lower-case
%    Gamma.  Gamma is defined as
%
%      Gamma = limit ( M -> Infinity )
%        ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
%
%  Modified:
%
%    08 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Output, real VALUE, the value of the
%    Euler-Mascheroni constant.
%
  value = 0.5772156649015328;
