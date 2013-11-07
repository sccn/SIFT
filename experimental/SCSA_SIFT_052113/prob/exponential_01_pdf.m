function pdf = exponential_01_pdf ( x )

%% EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF.
%
%  Formula:
%
%    PDF(X) = EXP ( - X )
%
%  Modified:
%
%    09 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%    0.0D+00 <= X
%
%    Output, real PDF, the value of the PDF.
%
  if ( x < 0.0 )
    pdf = 0.0;
  else
    pdf = exp ( - x );
  end
