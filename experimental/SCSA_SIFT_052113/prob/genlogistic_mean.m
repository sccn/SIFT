function mean = genlogistic_mean ( a, b, c )

%% GENLOGISTIC_MEAN returns the mean of the Generalized Logistic PDF.
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
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B,
%    0.0 < C.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = a + b * ( euler_constant ( 'DUMMY' ) + digamma ( c ) );
