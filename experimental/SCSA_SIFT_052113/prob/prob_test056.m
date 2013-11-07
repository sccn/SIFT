%% TEST056 tests EMPIRICAL_DISCRETE_CDF.
%% TEST056 tests EMPIRICAL_DISCRETE_PDF.
%
  clear

  a = 6;

  b(1:a) = [ 1.0, 1.0, 3.0, 2.0, 1.0,  2.0 ];
  c(1:a) = [ 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 ];

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST056\n' );
  fprintf ( 1, '  For the Empirical Discrete PDF.\n' );
  fprintf ( 1, '  EMPIRICAL_DISCRETE_PDF evaluates the PDF.\n' );
  fprintf ( 1, '  EMPIRICAL_DISCRETE_CDF evaluates the CDF.\n' );

  check = empirical_discrete_check ( a, b, c );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST056 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A = %6d\n', a );
  r8vec_print ( a, b, '  PDF parameter B:' );
  r8vec_print ( a, c, '  PDF parameter C:' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '      X      PDF(X)      CDF(X)\n' );
  fprintf ( 1, '\n' );

  for i = -2 : 12
    x = i;
    pdf = empirical_discrete_pdf ( x, a, b, c );
    cdf = empirical_discrete_cdf ( x, a, b, c );
    fprintf ( 1, '  %8f  %14f  %14f\n', x, pdf, cdf );
  end
