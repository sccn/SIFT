function x = cosine_cdf_inv ( cdf, a, b )

%% COSINE_CDF_INV inverts the Cosine CDF.
%
%  Discussion:
%
%    A simple bisection method is used on the interval
%    [ A - PI * B, A + PI * B ].
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
%    Input, real CDF, the value of the CDF.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real X, the corresponding argument of the CDF.
%
  it_max = 100;
  tol = 0.0001;

  if ( cdf <= 0.0 )
    x = a - pi * b;
    return
  elseif ( 1.0 <= cdf )
    x = a + pi * b;
    return
  end

  x1 = a - pi * b;
  cdf1 = 0.0;

  x2 = a + pi * b;
  cdf2 = 1.0;
%
%  Now use bisection.
%
  it = 0;

  for it = 1 : it_max

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = cosine_cdf ( x3, a, b );

    if ( abs ( cdf3 - cdf ) < tol )
      x = x3;
      return
    end

    if ( ( cdf3 < cdf & cdf1 < cdf ) | ( cdf < cdf3 & cdf < cdf1 ) )
      x1 = x3;
      cdf1 = cdf3;
    else
      x2 = x3;
      cdf2 = cdf3;
    end

  end

  fprintf ( 1, '\n' );
  fprintf ( 1, 'COSINE_CDF_INV - Fatal error!\n' );
  fprintf ( 1, '  Iteration limit exceeded.\n' );
