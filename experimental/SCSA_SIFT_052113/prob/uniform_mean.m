function mean = uniform_mean ( a, b )

%% UNIFORM_MEAN returns the mean of the Uniform PDF.
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
%    Input, real A, B, the parameters of the PDF.
%    A < B.
%
%    Output, real MEAN, the mean of the discrete uniform PDF.
%
  mean = 0.5 * ( a + b );
