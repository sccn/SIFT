function cdf = gamma_cdf ( x, a, b, c )

%% GAMMA_CDF evaluates the Gamma CDF.
%
%  Modified:
%
%    11 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%    A <= X
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real CDF, the value of the CDF.
%
  x2 = ( x - a ) / b;
  p2 = c;

  cdf = gamma_inc ( p2, x2 );
