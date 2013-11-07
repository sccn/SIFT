%% TEST104 tests LOG_UNIFORM_CDF;
%% TEST104 tests LOG_UNIFORM_INV.
%% TEST104 tests LOG_UNIFORM_PDF;
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST104\n' );
  fprintf ( 1, '  For the Log Uniform PDF:\n' );
  fprintf ( 1, '  LOG_UNIFORM_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  LOG_UNIFORM_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  LOG_UNIFORM_PDF evaluates the PDF;\n' );

  a = 2.0;
  b = 20.0;

  check = log_uniform_check ( a, b );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST104 - Fatal error!\n' );
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

    [ x, seed ] = log_uniform_sample ( a, b, seed );

    pdf = log_uniform_pdf ( x, a, b );

    cdf = log_uniform_cdf ( x, a, b );

    x2 = log_uniform_cdf_inv ( cdf, a, b );

    fprintf ( 1, ' %14f  %14f  %14f  %14f\n', x, pdf, cdf, x2 );

  end
