function pdf = circular_normal_01_pdf ( x, pdf )

%% CIRCULAR_NORMAL_01_PDF evaluates the Circular Normal 01 PDF.
%
%  Formula:
%
%    PDF(X) = EXP ( - 0.5D+00 * ( X(1)**2 + X(2)**2 ) ) / ( 2 * PI )
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
%    Input, real X(2), the argument of the PDF.
%
%    Output, real PDF, the value of the PDF.
%
  pdf = exp ( - 0.5 * ( x(1)^2 + x(2)^2 ) ) / ( 2.0 * pi );
