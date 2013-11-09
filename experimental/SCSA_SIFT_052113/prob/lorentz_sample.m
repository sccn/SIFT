function [ x, seed ] = lorentz_sample ( seed )

%% LORENTZ_SAMPLE samples the Lorentz PDF.
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
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = lorentz_cdf_inv ( cdf );
