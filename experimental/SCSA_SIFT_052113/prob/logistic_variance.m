function variance = logistic_variance ( a, b )

%% LOGISTIC_VARIANCE returns the variance of the Logistic PDF.
%
%  Modified:
%
%    13 September 2004
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
  variance = ( pi * b )^2 / 3.0;
