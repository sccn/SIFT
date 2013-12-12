%% TEST1555 tests VON_MISES_CDF.
%% TEST1555 tests VON_MISES_CDF_VALUES.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST1555:\n' );
  fprintf ( 1, '  VON_MISES_CDF evaluates the von Mises CDF.\n' );
  fprintf ( 1, '  VON_MISES_CDF_VALUES returns some exact values.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  A is the dominant angle;\n' );
  fprintf ( 1, '  B is a measure of spread;\n' );
  fprintf ( 1, '  X is the angle;\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '      A     B         X   Exact F     Computed F\n' );
  fprintf ( 1, '\n' );
  n_data = 0;

  while ( 1 )

    [ n_data, a, b, x, fx ] = von_mises_cdf_values ( n_data );

    if ( n_data == 0 );
      break
    end

    fx2 = von_mises_cdf ( x, a, b );

    fprintf ( 1, '  %8f  %8f  %8f  %14f  %14f\n', a, b, x, fx, fx2 );

  end
