function value = r8_ceiling ( x )

%% R8_CEILING returns the "ceiling" of a number.
%
%  Definition:
%
%    The "ceiling" of X is the value of X rounded towards plus infinity.
%
%  Examples:
%
%    X         Value
%
%   -1.1      -1.0
%   -1.0      -1.0
%   -0.9       0.0
%   -0.1       0.0
%    0.0       0.0
%    0.1       1.0
%    0.9       1.0
%    1.0       1.0
%    1.1       2.0
%    2.9       3.0
%    3.0       3.0
%    3.14159   4.0
%
%  Modified:
%
%    25 May 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the number whose ceiling is desired.
%
%    Output, real VALUE, the ceiling of X.
%
  value = ceil ( x );
