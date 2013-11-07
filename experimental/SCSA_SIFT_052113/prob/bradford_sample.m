function [ x, seed ] = bradford_sample ( a, b, c, seed )

%% BRADFORD_SAMPLE samples the Bradford PDF.
%
%  Modified:
%
%    03 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, C, the parameters of the PDF.
%    A < B,
%    0.0 < C.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%
%    Output, integer SEED, a seed for the random number generator.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = a + ( b - a ) * ( ( c + 1.0 )^cdf - 1.0 ) / c;
