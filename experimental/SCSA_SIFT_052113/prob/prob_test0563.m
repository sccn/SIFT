%% TEST0563 tests ENGLISH_SENTENCE_LENGTH_CDF, ENGLISH_SENTENCE_LENGTH_CDF_INV and ENGLISH_SENTENCE_LENGTH_PDF.
%
%  Modified:
%
%    02 August 2006
%
%  Author:
%
%    John Burkardt
%
  clear

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST0565\n' );
  fprintf ( 1, '  For the English Sentence Length PDF:\n' );
  fprintf ( 1, '  ENGLISH_SENTENCE_LENGTH_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  ENGLISH_SENTENCE_LENGTH_PDF evaluates the PDF;\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = english_sentence_length_sample ( seed );

    pdf = english_sentence_length_pdf ( x );

    cdf = english_sentence_length_cdf ( x );

    x2 = english_sentence_length_cdf_inv ( cdf );

    fprintf ( 1, '  %12d  %12f  %12f  %12d\n', x, pdf, cdf, x2 );

  end
