function pdf = gumbel_pdf ( x )

%% GUMBEL_PDF evaluates the Gumbel PDF.
%
%  Formula:
%
%    PDF(X) = EXP ( - X - EXP ( - X  ) ).
%
%  Discussion:
%
%    GUMBEL_PDF(X) = EXTREME_PDF(X)(0,1)
%
%  Modified:
%
%    12 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Eric Weisstein, editor,
%    CRC Concise Encylopedia of Mathematics,
%    CRC Press, 1998.
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%
%    Output, real PDF, the value of the PDF.
%
  pdf = exp ( - x - exp ( - x ) );
