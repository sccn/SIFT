function pdf = weibull_discrete_pdf ( x, a, b )

%% WEIBULL_DISCRETE_PDF evaluates the discrete Weibull PDF.
%
%  Formula:
%
%    PDF(X)(A,B) = ( 1 - A )**X**B - ( 1 - A )**(X+1)**B.
%
%  Discussion:
%
%    WEIBULL_DISCRETE_PDF(X)(A,1) = GEOMETRIC_PDF(X)(A)
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
%    Input, integer X, the argument of the PDF.
%    0 <= X
%
%    Input, real A, B, the parameters that define the PDF.
%    0 <= A <= 1,
%    0 < B.
%
%    Output, real PDF, the value of the PDF.
%
  if ( x < 0 )
    pdf = 0.0;
  else
    pdf = ( 1.0 - a )^(x^b) - ( 1.0 - a )^((x+1)^b);
  end
