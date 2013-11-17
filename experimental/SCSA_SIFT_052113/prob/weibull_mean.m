function mean = weibull_mean ( a, b, c )

%% WEIBULL_MEAN returns the mean of the Weibull PDF.
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
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = b * gamma ( ( c + 1.0 ) / c ) + a;
