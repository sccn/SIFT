function [ x, seed ] = deranged_sample ( a, seed )

%% DERANGED_SAMPLE samples the Deranged PDF.
%
%  Modified:
%
%    07 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer A, the number of items.
%    1 <= A.
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, integer X, a sample of the PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = deranged_cdf_inv ( cdf, a );
