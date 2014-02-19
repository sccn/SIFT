function mean = uniform_discrete_mean ( a, b )

%% UNIFORM_DISCRETE_MEAN returns the mean of the Uniform discrete PDF.
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
%    Input, integer A, B, the parameters of the PDF.
%    A <= B.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = 0.5 * ( a + b );
