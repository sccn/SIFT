function mean = runs_mean ( m, n )

%% RUNS_MEAN returns the mean of the Runs PDF.
%
%  Modified:
%
%    23 January 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the parameters of the PDF.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = ( m + 2 * m * n + n ) / ( m + n );

