%% TEST115 tests NORMAL_01_CDF;
%% TEST115 tests NORMAL_01_CDF_INV.
%% TEST115 tests NORMAL_01_PDF;
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST115\n' );
  fprintf ( 1, '  For the Normal 01 PDF:\n' );
  fprintf ( 1, '  NORMAL_01_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  NORMAL_01_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  NORMAL_01_PDF evaluates the PDF;\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  seed = 123456789;

  for i = 1 : 10

    [ x, seed ] = normal_01_sample ( seed );

    pdf = normal_01_pdf ( x );

    cdf = normal_01_cdf ( x );

    x2 = normal_01_cdf_inv ( cdf );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
