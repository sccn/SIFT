function [ x, seed ] = uniform_discrete_sample ( a, b, seed )

%% UNIFORM_DISCRETE_SAMPLE samples the Uniform discrete PDF.
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
%    Input, integer A, B, the parameters of the PDF.
%    A <= B.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, integer X, a sample of the PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = uniform_discrete_cdf_inv ( cdf, a, b );
