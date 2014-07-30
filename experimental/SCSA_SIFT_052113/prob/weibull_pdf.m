function pdf = weibull_pdf ( x, a, b, c )

%% WEIBULL_PDF evaluates the Weibull PDF.
%
%  Formula:
%
%    PDF(X)(A,B,C) = ( C / B ) * ( ( X - A ) / B )**( C - 1 )
%     * EXP ( - ( ( X - A ) / B )**C ).
%
%  Discussion:
%
%    The Weibull PDF is also known as the Frechet PDF.
%
%    WEIBULL_PDF(X)(A,B,1) is the Exponential PDF.
%
%    WEIBULL_PDF(X)(0,1,2) is the Rayleigh PDF.
%
%  Modified:
%
%    22 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%    A <= X
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real PDF, the value of the PDF.
%
  if ( x < a )

    pdf = 0.0;

  else

    y = ( x - a ) / b;

    pdf = ( c / b ) * y^( c - 1.0 )  / exp ( y^c );

  end
