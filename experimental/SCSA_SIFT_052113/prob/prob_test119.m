%% TEST119 tests PARETO_CDF.
%% TEST119 tests PARETO_CDF_INV.
%% TEST119 tests PARETO_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST119\n' );
  fprintf ( 1, '  For the Pareto PDF:\n' );
  fprintf ( 1, '  PARETO_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  PARETO_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  PARETO_PDF evaluates the PDF;\n' );

  a = 2.0;
  b = 3.0;

  check = pareto_check ( a, b );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST119 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A =             %14f\n', a );
  fprintf ( 1, '  PDF parameter B =             %14f\n', b );

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = pareto_sample ( a, b, seed );

    pdf = pareto_pdf ( x, a, b );

    cdf = pareto_cdf ( x, a, b );

    x2 = pareto_cdf_inv ( cdf, a, b );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
