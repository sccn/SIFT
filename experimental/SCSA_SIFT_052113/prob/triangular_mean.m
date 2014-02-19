function mean = triangular_mean ( a, b )

%% TRIANGULAR_MEAN returns the mean of the Triangular PDF.
%
%  Modified:
%
%    13 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    A < B.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = 0.5 * ( a + b );

