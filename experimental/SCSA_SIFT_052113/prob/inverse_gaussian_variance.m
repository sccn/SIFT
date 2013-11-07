function variance = inverse_gaussian_variance ( a, b )

%% INVERSE_GAUSSIAN_VARIANCE returns the variance of the Inverse Gaussian PDF.
%
%  Modified:
%
%    12 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < A,
%    0.0 < B.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = a^3 / b;
