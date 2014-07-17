function x = sech_cdf_inv ( cdf, a, b )

%% SECH_CDF_INV inverts the Hyperbolic Secant CDF.
%
%  Modified:
%
%    20 September 2004
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
  if ( cdf <= 0.0 )
    x = - r8_huge ( x );
  elseif ( cdf < 1.0 )
    x = a + b * log ( tan ( 0.5 * pi * cdf ) );
  elseif ( 1.0 <= cdf )
    x = r8_huge ( x );
  end
