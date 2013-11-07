%% TEST059 tests ERROR_F.
%
%  Modified:
%
%    17 November 2006
%
%  Author:
%
%    John Burkardt
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST059\n' );
  fprintf ( 1, '  ERROR_F evaluates the error function ERF.\n' );

  x = 1.0;

  erfx = error_f ( x );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  ERF argument X = %14f\n', x );
  fprintf ( 1, '  ERF value        %14f\n', erfx );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  (Expected answer is 0.843)\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Test:\n' );
  fprintf ( 1, '    0.5 * ( ERF(X/SQRT(2)) + 1 ) = Normal_CDF(X)\n' );
  fprintf ( 1, '\n' );

  x = 1.0;
  x2 = x / sqrt ( 2.0 );
  erfx = error_f ( x2 );

  cdf = normal_01_cdf ( x );

  fprintf ( 1, '  0.5 * ( ERF(X/SQRT(2)) + 1 ) = %14f\n', ...
    0.5 * ( erfx + 1.0D+00 ) );
  fprintf ( 1, '  Normal_CDF(X) = %14f\n', cdf );
