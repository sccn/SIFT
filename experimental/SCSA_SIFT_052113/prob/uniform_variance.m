function variance = uniform_variance ( a, b )

%% UNIFORM_VARIANCE returns the variance of the Uniform PDF.
%
%  Modified:
%
%    22 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    A < B.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = ( b - a )^2 / 12.0;
