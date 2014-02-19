function pdf = cauchy_pdf ( x, a, b )

%% CAUCHY_PDF evaluates the Cauchy PDF.
%
%  Formula:
%
%    PDF(X)(A,B) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )**2 ) )
%
%  Discussion:
%
%    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
%    has some unusual properties.  In particular, the integrals for the
%    expected value and higher order moments are "singular", in the
%    sense that the limiting values do not exist.  A result can be
%    obtained if the upper and lower limits of integration are set
%    equal to +T and -T, and the limit as T=>INFINITY is taken, but
%    this is a very weak and unreliable sort of limit.
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
%    Input, real X, the argument of the PDF.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real PDF, the value of the PDF.
%
  y = ( x - a ) / b;

  pdf = 1.0 / ( pi * b * ( 1.0 + y * y ) );
