function pdf = logistic_pdf ( x, a, b )

%% LOGISTIC_PDF evaluates the Logistic PDF.
%
%  Formula:
%
%    PDF(X)(A,B) = EXP ( ( A - X ) / B ) /
%      ( B * ( 1 + EXP ( ( A - X ) / B ) )**2 )
%
%  Discussion:
%
%    The Logistic PDF is also known as the Sech-Squared PDF.
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
  temp = exp ( ( a - x ) / b );

  pdf = temp / ( b * ( 1.0D+00 + temp )^2 );
