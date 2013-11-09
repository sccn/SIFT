function [ x, seed ] = english_word_length_sample ( seed )

%% ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Length PDF.
%
%  Modified:
%
%    29 July 2006
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Henry Kucera, Winthrop Francis,
%    Computational Analysis of Present-Day American English,
%    Brown University Press, 1967.
%
%  Parameters:
%
%    Input/output, integer SEED, a seed for the random number generator.
%
%    Output, integer X, a sample of the PDF.
%
  [ cdf, seed ] = r8_uniform_01 ( seed );

  x = english_word_length_cdf_inv ( cdf );
