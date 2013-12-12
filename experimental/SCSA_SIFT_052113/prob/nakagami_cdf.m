function cdf = nakagami_cdf ( x, a, b, c )

%% NAKAGAMI_CDF evaluates the Nakagami CDF.
%
%  Modified:
%
%    17 September 2004
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
%    0.0 < B
%    0.0 < C.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= 0.0 )

    cdf = 0.0;

  elseif ( 0.0 < x )

    y = ( x - a ) / b;
    x2 = c * y * y;
    p2 = c;

    cdf = gamma_inc ( p2, x2 );

  end
