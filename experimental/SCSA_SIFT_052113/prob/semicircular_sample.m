function [ x, seed ] = semicircular_sample ( a, b, seed )

%% SEMICIRCULAR_SAMPLE samples the Semicircular PDF.
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
%    Input, real A, B, the parameters of the PDF.
%    0.0 < B.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%
%    Output, integer SEED, a seed for the random number generator.
%
  [ radius, seed ] = r8_uniform_01 ( seed );
  radius = b * sqrt ( radius );
  [ angle, seed ] = r8_uniform_01 ( seed );
  x = a + radius * cos ( pi * angle );
