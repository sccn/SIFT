function variance = lorentz_variance ( dummy )

%% LORENTZ_VARIANCE returns the variance of the Lorentz PDF.
%
%  Discussion:
%
%    The variance of the Lorentz PDF is not well defined.  This routine
%    is made available for completeness only, and simply returns
%    a "very large" number.
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
%    Output, real VARIANCE, the mean of the PDF.
%
  variance = r8_huge ( 'DUMMY' );
