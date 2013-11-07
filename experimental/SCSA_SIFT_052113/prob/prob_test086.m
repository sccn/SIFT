%% TEST086 tests GUMBEL_CDF;
%% TEST086 tests GUMBEL_CDF_INV.
%% TEST086 tests GUMBEL_PDF;
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST086\n' );
  fprintf ( 1, '  For the Gumbel PDF:\n' );
  fprintf ( 1, '  GUMBEL_CDF evaluates the CDF.\n' );
  fprintf ( 1, '  GUMBEL_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  GUMBEL_PDF evaluates the PDF.\n' );

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = gumbel_sample ( seed );

    pdf = gumbel_pdf ( x );

    cdf = gumbel_cdf ( x );

    x2 = gumbel_cdf_inv ( cdf );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
