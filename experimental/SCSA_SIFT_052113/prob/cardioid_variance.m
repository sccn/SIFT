function variance = cardioid_variance ( x, a, b )

%% CARDIOID_VARIANCE returns the variance of the Cardioid PDF.
%
%  Modified:
%
%    31 July 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 <= B <= 0.5.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = a;
