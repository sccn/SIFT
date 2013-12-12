function mean = chi_mean ( a, b, c )

%% CHI_MEAN returns the mean of the Chi PDF.
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
%    Output, real MEAN, the mean value.
%
  mean = a + sqrt ( 2.0 ) * b * gamma ( 0.5 * ( c + 1.0 ) ) / gamma ( 0.5 * c );
