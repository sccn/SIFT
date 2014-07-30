function pdf = genlogistic_pdf ( x, a, b, c )

%% GENLOGISTIC_PDF evaluates the Generalized Logistic PDF.
%
%  Formula:
%
%    PDF(X)(A,B,C) = ( C / B ) * EXP ( ( A - X ) / B ) /
%      ( ( 1 + EXP ( ( A - X ) / B ) )**(C+1) )
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
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real PDF, the value of the PDF.
%
  y = ( x - a ) / b;

  pdf = ( c / b ) * exp ( - y ) / ( 1.0 + exp ( - y ) )^( c + 1.0 );
