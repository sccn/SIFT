function pdf = zipf_pdf ( x, a )

%% ZIPF_PDF evaluates the Zipf PDF.
%
%  Formula:
%
%    PDF(X)(A) = ( 1 / X**A ) / C
%
%    where the normalizing constant is chosen so that
%
%    C = Sum ( 1 <= I < Infinity ) 1 / I**A.
%
%  Discussion:
%
%    From observation, the frequency of different words in long
%    sequences of text seems to follow the Zipf PDF, with
%    parameter A slightly greater than 1.  The Zipf PDF is sometimes
%    known as the "discrete Pareto" PDF.
%
%  Modified:
%
%    08 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X, the argument of the PDF.
%    1 <= N
%
%    Input, real A, the parameter of the PDF.
%    1.0 < A.
%
%    Output, real PDF, the value of the PDF.
%
  if ( x < 1 )

    pdf = 0.0;

  else

    c = zeta ( a );
    pdf = ( 1.0 / x^a ) / c;

  end
