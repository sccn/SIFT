%% TEST033 tests CHI_SQUARE_CDF.
%% TEST033 tests CHI_SQUARE_CDF_INV.
%% TEST033 tests CHI_SQUARE_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST033\n' );
  fprintf ( 1, '  For the central chi square PDF:\n' );
  fprintf ( 1, '  CHI_SQUARE_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  CHI_SQUARE_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  CHI_SQUARE_PDF evaluates the PDF;\n' );

  a = 4.0;

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A = %14f\n', a );

  check = chi_square_check ( a );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST033 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = chi_square_sample ( a, seed );

    pdf = chi_square_pdf ( x, a );

    cdf = chi_square_cdf ( x, a );

    x2 = chi_square_cdf_inv ( cdf, a );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
