%% TEST073 tests FOLDED_NORMAL_CDF.
%% TEST073 tests FOLDED_NORMAL_CDF_INV.
%% TEST073 tests FOLDED_NORMAL_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST073\n' );
  fprintf ( 1, '  For the Folded Normal PDF:\n' );
  fprintf ( 1, '  FOLDED_NORMAL_CDF evaluates the CDF.\n' );
  fprintf ( 1, '  FOLDED_NORMAL_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  FOLDED_NORMAL_PDF evaluates the PDF.\n' );

  a = 2.0;
  b = 3.0;

  check = folded_normal_check ( a, b );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST073 - Fatal error!\n' );
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

    [ x, seed ] = folded_normal_sample ( a, b, seed );

    pdf = folded_normal_pdf ( x, a, b );

    cdf = folded_normal_cdf ( x, a, b );

    x2 = folded_normal_cdf_inv ( cdf, a, b );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
