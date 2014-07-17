function mean = exponential_mean ( a, b )

%% EXPONENTIAL_MEAN returns the mean of the Exponential PDF.
%
%  Modified:
%
%    09 September 2004
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
%    Output, real MEAN, the mean of the PDF.
%
  mean = a + b;

