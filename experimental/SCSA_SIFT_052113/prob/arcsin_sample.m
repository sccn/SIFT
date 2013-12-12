function [ x, seed ] = arcsin_sample ( a, seed )

%% ARCSIN_SAMPLE samples the Arcsin PDF.
%
%  Modified:
%
%    01 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, the parameter of the CDF.
%    A must be positive.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = arcsin_cdf_inv ( cdf, a );
