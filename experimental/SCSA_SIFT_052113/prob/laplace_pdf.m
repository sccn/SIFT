function pdf = laplace_pdf ( x, a, b )

%% LAPLACE_PDF evaluates the Laplace PDF.
%
%  Formula:
%
%    PDF(X)(A,B) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
%
%  Discussion:
%
%    The Laplace PDF is also known as the Double Exponential PDF.
%
%  Modified:
%
%    13 September 2004
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
  pdf = exp ( - abs ( x - a ) / b ) / ( 2.0 * b );
