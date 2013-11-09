function mean = maxwell_mean ( a )

%% MAXWELL_MEAN returns the mean of the Maxwell PDF.
%
%  Modified:
%
%    16 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, the parameter of the PDF.
%    0 < A.
%
%    Output, real MEAN, the mean value.
%
  mean = sqrt ( 2.0 ) * a * gamma ( 2.0 ) / gamma ( 1.5 );
