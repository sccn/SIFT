function variance = triangular_variance ( a, b )

%% TRIANGULAR_VARIANCE returns the variance of the Triangular PDF.
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
  variance = ( b - a )^2 / 24.0;
