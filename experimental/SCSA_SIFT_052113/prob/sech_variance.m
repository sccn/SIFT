function variance = sech_variance ( a, b )

%% SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
%
%  Modified:
%
%    20 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = 0.25 * pi * pi * b * b;
