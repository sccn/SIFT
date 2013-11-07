function x = cauchy_cdf_inv ( cdf, a, b )

%% CAUCHY_CDF_INV inverts the Cauchy CDF.
%
%  Modified:
%
%    05 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real CDF, the value of the CDF.
%    0.0 <= CDF <= 1.0.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real X, the corresponding argument.
%
  x = a + b * tan ( pi * ( cdf - 0.5 ) );
