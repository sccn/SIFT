function cdf = chi_cdf ( x, a, b, c )

%% CHI_CDF evaluates the Chi CDF.
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
%    Input, real X, the argument of the PDF.
%
%    Input, real A, B, C, the parameters of the PDF.
%    0 < B,
%    0 < C.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= a )

    cdf = 0.0;

  else

    y = ( x - a ) / b;
    x2 = 0.5 * y * y;
    p2 = 0.5 * c;

    cdf = gamma_inc ( p2, x2 );

  end
