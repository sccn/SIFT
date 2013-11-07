function mean = inverse_gaussian_mean ( a, b )

%% INVERSE_GAUSSIAN_MEAN returns the mean of the Inverse Gaussian PDF.
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
%    Output, real MEAN, the mean of the PDF.
%
  mean = a;
