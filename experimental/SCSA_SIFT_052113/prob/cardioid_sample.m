function [ x, seed ] = cardioid_sample ( a, b, seed )

%% CARDIOID_SAMPLE samples the Cardioid PDF.
%
%  Modified:
%
%    30 July 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    0.0 <= B <= 0.5.
%
%    Input/output, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the PDF.
%    A - PI <= X <= A + PI.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = cardioid_cdf_inv ( cdf, a, b );
