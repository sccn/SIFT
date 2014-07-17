function [ x, seed ] = f_sample ( m, n, seed )

%% F_SAMPLE samples the F central PDF.
%
%  Modified:
%
%    30 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer M, N, the parameters of the PDF.
%    1 <= M,
%    1 <= N.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ xm, seed ] = chi_square_sample ( m, seed );

  [ xn, seed ] = chi_square_sample ( n, seed );

  x = n * xm / ( m * xn );
