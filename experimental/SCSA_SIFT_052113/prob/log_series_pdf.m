function pdf = log_series_pdf ( x, a )

%% LOG_SERIES_PDF evaluates the Logarithmic Series PDF.
%
%  Formula:
%
%    PDF(X)(A) = - A**X / ( X * log ( 1 - A ) )
%
%  Modified:
%
%    20 February 1999
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X, the argument of the PDF.
%    0 < X
%
%    Input, real A, the parameter of the PDF.
%    0.0 < A < 1.0.
%
%    Output, real PDF, the value of the PDF.
%
  if ( x <= 0 )
    pdf = 0.0;
  else
    pdf = - a^x / ( x * log ( 1.0 - a ) );
  end
