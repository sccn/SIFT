function variance = negative_binomial_variance ( a, b )

%% NEGATIVE_BINOMIAL_VARIANCE returns the variance of the Negative Binomial PDF.
%
%  Modified:
%
%    19 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer A, a parameter of the PDF.
%    0 <= A.
%
%    Input, real B, a parameter of the PDF.
%    0 < B <= 1.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = a * ( 1.0 - b ) / ( b * b );
