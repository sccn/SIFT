function mean = reciprocal_mean ( a, b )

%% RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
%
%  Modified:
%
%    20 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < A <= B.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = ( a - b ) / log ( a / b );
