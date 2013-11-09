function variance = chi_variance ( a, b, c )

%% CHI_VARIANCE returns the variance of the Chi PDF.
%
%  Modified:
%
%    06 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, C, the parameters of the PDF.
%    0 < B,
%    0 < C.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = b * b * ( c - 2.0 * ( gamma ( 0.5 * ( c + 1.0 ) ) / gamma ( 0.5 * c ) )^2 );
