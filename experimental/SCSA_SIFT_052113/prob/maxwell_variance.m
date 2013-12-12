function variance = maxwell_variance ( a )

%% MAXWELL_VARIANCE returns the variance of the Maxwell PDF.
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
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = a * a * ( 3.0 - 2.0 * ( gamma ( 2.0 ) / gamma ( 1.5 ) )^2 );
