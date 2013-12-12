function [ x, seed ] = english_sentence_length_sample ( seed )

%% ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length PDF.
%
%  Modified:
%
%    02 August 2006
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

  x = english_sentence_length_cdf_inv ( cdf );
