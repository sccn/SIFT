function mean = triangle_mean ( a, b, c )

%% TRIANGLE_MEAN returns the mean of the Triangle PDF.
%
%  Modified:
%
%    21 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, C, the parameters of the PDF.
%    A <= B <= C and A < C.
%
%    Output, real MEAN, the mean of the discrete uniform PDF.
%
  mean = a + ( c + b - 2.0 * a ) / 3.0;
