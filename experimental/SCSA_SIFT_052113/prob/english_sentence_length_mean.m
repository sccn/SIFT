function mean = english_sentence_length_mean ( DUMMY )

%% ENGLISH_SENTENCE_LENGTH_MEAN evaluates the mean of the English Sentence Length PDF.
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
%    Output, real MEAN, the mean of the PDF.
%
  word_length_max = 79;

  pdf_vec = [ ...
    0.00806, ...
    0.01370, ...
    0.01862, ...
    0.02547, ...
    0.03043, ...
    0.03189, ...
    0.03516, ...
    0.03545, ...
    0.03286, ...
    0.03533, ...
    0.03562, ...
    0.03788, ...
    0.03669, ...
    0.03751, ...
    0.03518, ...
    0.03541, ...
    0.03434, ...
    0.03305, ...
    0.03329, ...
    0.03103, ...
    0.02867, ...
    0.02724, ...
    0.02647, ...
    0.02526, ...
    0.02086, ...
    0.02178, ...
    0.02128, ...
    0.01801, ...
    0.01690, ...
    0.01556, ...
    0.01512, ...
    0.01326, ...
    0.01277, ...
    0.01062, ...
    0.01051, ...
    0.00901, ...
    0.00838, ...
    0.00764, ...
    0.00683, ...
    0.00589, ...
    0.00624, ...
    0.00488, ...
    0.00477, ...
    0.00406, ...
    0.00390, ...
    0.00350, ...
    0.00318, ...
    0.00241, ...
    0.00224, ...
    0.00220, ...
    0.00262, ...
    0.00207, ...
    0.00174, ...
    0.00174, ...
    0.00128, ...
    0.00121, ...
    0.00103, ...
    0.00117, ...
    0.00124, ...
    0.00082, ...
    0.00088, ...
    0.00061, ...
    0.00061, ...
    0.00075, ...
    0.00063, ...
    0.00056, ...
    0.00052, ...
    0.00057, ...
    0.00031, ...
    0.00029, ...
    0.00021, ...
    0.00017, ...
    0.00021, ...
    0.00034, ...
    0.00031, ...
    0.00011, ...
    0.00011, ...
    0.00008, ...
    0.00006 ];
  pdf_sum = 0.99768;

  mean = 0.0;
  for j = 1 : word_length_max
    mean = mean + j * pdf_vec(j);
  end

  mean = mean / pdf_sum;
