function mean = extreme_values_mean ( a, b )

%% EXTREME_VALUES_MEAN returns the mean of the Extreme Values PDF.
%
%  Modified:
%
%    09 September 2004
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
  mean = a + b * euler_constant ( 'DUMMY' );
