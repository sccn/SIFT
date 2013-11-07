function variance = log_series_variance ( a )

%% LOG_SERIES_VARIANCE returns the variance of the Logarithmic Series PDF.
%
%  Modified:
%
%    20 February 1999
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, the parameter of the PDF.
%    0.0 < A < 1.0.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  alpha = - 1.0 / log ( 1.0 - a );

  variance = a * alpha * ( 1.0 - alpha * a ) / ( 1.0 - a )^2;

