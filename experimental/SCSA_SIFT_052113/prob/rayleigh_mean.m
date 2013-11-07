function mean = rayleigh_mean ( a )

%% RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
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
%    Input, real A, the parameter of the PDF.
%    0.0 < A.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = a * sqrt ( 0.5 * pi );
