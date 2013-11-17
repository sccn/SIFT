function variance = half_normal_variance ( a, b, variance )

%% HALF_NORMAL_VARIANCE returns the variance of the Half Normal PDF.
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
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = b * b * ( 1.0 - 2.0 / pi );
