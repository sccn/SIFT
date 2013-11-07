function variance = reciprocal_variance ( a, b )

%% RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
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
%    Input, real A, B, the parameters of the PDF.
%    0.0 < A <= B.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  d = log ( a / b );

  variance = ( a - b )* ( a * ( d - 2.0 ) + b * ( d + 2.0 ) ) / ( 2.0 * d * d );
