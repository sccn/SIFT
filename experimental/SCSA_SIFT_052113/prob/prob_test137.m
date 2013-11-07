%% TEST137 tests SEMICIRCULAR_CDF.
%% TEST137 tests SEMICIRCULAR_CDF_INV.
%% TEST137 tests SEMICIRCULAR_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST137\n' );
  fprintf ( 1, '  For the Semicircular PDF:\n' );
  fprintf ( 1, '  SEMICIRCULAR_CDF evaluates the CDF.\n' );
  fprintf ( 1, '  SEMICIRCULAR_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  SEMICIRCULAR_PDF evaluates the PDF.\n' );

  a = 3.0;
  b = 2.0;

  check = semicircular_check ( a, b );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST137 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A =         %14f\n', a );
  fprintf ( 1, '  PDF parameter B =         %14f\n', b );

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = semicircular_sample ( a, b, seed );

    pdf = semicircular_pdf ( x, a, b );

    cdf = semicircular_cdf ( x, a, b );

    x2 = semicircular_cdf_inv ( cdf, a, b );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
