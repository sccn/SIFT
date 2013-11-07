function variance = log_normal_variance ( a, b )

%% LOG_NORMAL_VARIANCE returns the variance of the Lognormal PDF.
%
%  Modified:
%
%    12 February 1999
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
  variance = exp ( 2.0 * a + b * b ) * ( exp ( b * b ) - 1.0 );
