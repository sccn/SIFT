function mean = log_series_mean ( a )

%% LOG_SERIES_MEAN returns the mean of the Logarithmic Series PDF.
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
%    Output, real MEAN, the mean of the PDF.
%
  mean = - a / ( ( 1.0 - a ) * log ( 1.0 - a ) );
