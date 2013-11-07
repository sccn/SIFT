function pdf = gamma_pdf ( x, a, b, c )

%% GAMMA_PDF evaluates the Gamma PDF.
%
%  Formula:
%
%    PDF(X)(A,B,C) = EXP ( - ( X - A ) / B ) * ( ( X - A ) / B )**(C-1)
%      / ( B * GAMMA ( C ) )
%
%  Discussion:
%
%    GAMMA_PDF(A,B,C), where C is an integer, is the Erlang PDF.
%    GAMMA_PDF(A,B,1) is the Exponential PDF.
%    GAMMA_PDF(0,2,C/2) is the Chi Squared PDF with C degrees of freedom.
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
%    A <= X.
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real PDF, the value of the PDF.
%
  if ( x <= a )

    pdf = 0.0;

  else

    y = ( x - a ) / b;

    pdf = y^( c - 1.0 ) / ( b * gamma ( c ) * exp ( y ) );

  end
