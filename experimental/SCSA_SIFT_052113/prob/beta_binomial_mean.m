function mean = beta_binomial_mean ( a, b, c )

%% BETA_BINOMIAL_MEAN returns the mean of the Beta Binomial PDF.
%
%  Modified:
%
%    02 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, parameters of the PDF.
%    0.0D+00 < A,
%    0.0D+00 < B.
%
%    Input, integer C, a parameter of the PDF.
%    0 <= N.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = c * a / ( a + b );

