%% TEST044 tests DIGAMMA.
%% TEST044 tests PSI_VALUES.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST044:\n' );
  fprintf ( 1, '  DIGAMMA evaluates the DIGAMMA or PSI function.\n' );
  fprintf ( 1, '  PSI_VALUES returns some exact values.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '     X       Exact F       DIGAMMA(X)\n' );
  fprintf ( 1, '\n' );

  n_data = 0;

  while ( 1 )

    [ n_data, x, fx ] = psi_values ( n_data );

    if ( n_data == 0 );
      break
    end

    if ( x <= 0.0 );
      continue;
    end

    fx2 = digamma ( x );

    fprintf ( 1, '  %8f  %14f  %14f\n', x, fx, fx2 );

  end
