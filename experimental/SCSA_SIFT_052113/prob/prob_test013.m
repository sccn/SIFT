%% TEST013 tests BETA_INC.
%% TEST013 tests BETA_INC_VALUES.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST013:\n' );
  fprintf ( 1, '  BETA_INC evaluates the normalized incomplete Beta\n' );
  fprintf ( 1, '    function BETA_INC(A,B,X).\n' );
  fprintf ( 1, '  BETA_INC_VALUES returns some exact values.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  A      B       X       Exact F       BETA_INC(A,B,X)\n' );
  fprintf ( 1, '\n' );
  n_data = 0;

  while ( 1 )

    [ n_data, a, b, x, fx ] = beta_inc_values ( n_data );

    if ( n_data == 0 );
      break
    end

    fx2 = beta_inc ( a, b, x );

    fprintf ( 1, '  %8f  %8f  %8f  %14f  %14f\n', a, b, x, fx, fx2 );

  end
