function pdf = normal_pdf ( x, a, b )

%% NORMAL_PDF evaluates the Normal PDF.
%
%  Formula:
%
%    PDF(X)(A,B)
%      = EXP ( - 0.5D+00 * ( ( X - A ) / B )**2 )
%      / ( B * SQRT ( 2 * PI ) )
%
%  Discussion:
%
%    The normal PDF is also known as the Gaussian PDF.
%
%  Modified:
%
%    17 September 2004
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

  pdf = exp ( - 0.5 * y * y )  / ( b * sqrt ( 2.0 * pi ) );
