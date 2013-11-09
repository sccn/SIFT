function mean = half_normal_mean ( a, b )

%% HALF_NORMAL_MEAN returns the mean of the Half Normal PDF.
%
%  Modified:
%
%    12 September 2004
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
%    Output, real MEAN, the mean of the PDF.
%
  mean = a + b * sqrt ( 2.0 / pi );
