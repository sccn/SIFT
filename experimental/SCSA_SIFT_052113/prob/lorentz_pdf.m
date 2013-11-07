function pdf = lorentz_pdf ( x )

%% LORENTZ_PDF evaluates the Lorentz PDF.
%
%  Formula:
%
%    PDF(X) = 1 / ( PI * ( 1 + X**2 ) )
%
%  Discussion:
%
%    The chief interest of the Lorentz PDF is that it is easily
%    inverted, and can be used to dominate other PDF's in an
%    acceptance/rejection method.
%
%    LORENTZ_PDF(X) = CAUCHY_PDF(X)(0,1)
%
%  Modified:
%
%    01 February 1999
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%
%    Output, real PDF, the value of the PDF.
%
  pdf = 1.0 / ( pi * ( 1.0 + x * x ) );
