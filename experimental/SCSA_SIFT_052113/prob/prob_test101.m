%% TEST101 tests LOG_SERIES_CDF.
%% TEST101 tests LOG_SERIES_CDF_INV.
%% TEST101 tests LOG_SERIES_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST101\n' );
  fprintf ( 1, '  For the Logseries PDF,\n' );
  fprintf ( 1, '  LOG_SERIES_CDF evaluates the CDF;\n' );
  fprintf ( 1, '  LOG_SERIES_CDF_INV inverts the CDF.\n' );
  fprintf ( 1, '  LOG_SERIES_PDF evaluates the PDF;\n' );

  a = 0.25;

  check = log_series_check ( a );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST101 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A =  %14f\n', a );

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, '       X            PDF           CDF            CDF_INV\n' );
  fprintf ( 1, '\n' );

  for i = 1 : 10

    [ x, seed ] = log_series_sample ( a, seed );

    pdf = log_series_pdf ( x, a );

    cdf = log_series_cdf ( x, a );

    x2 = log_series_cdf_inv ( cdf, a );

    fprintf ( 1, '  %14d  %14f  %14f  %14d\n', x, pdf, cdf, x2 );

  end
