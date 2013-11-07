function pdf = sech_pdf ( x, a, b )

%% SECH_PDF evaluates the Hypebolic Secant PDF.
%
%  Formula:
%
%    PDF(X)(A,B) = sech ( ( X - A ) / B ) / ( PI * B )
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
%    Input, real X, the argument of the PDF.
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Output, real PDF, the value of the PDF.
%
  y = ( x - a ) / b;

  pdf = sech ( y ) / ( pi * b );
