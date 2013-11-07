%% TEST1425 tests TFN.
%% TEST1425 tests OWEN_VALUES.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST1425\n' );
  fprintf ( 1, '  TFN evaluates Owen''s T function.\n' );
  fprintf ( 1, '  OWEN_VALUES returns some exact values.\n' );
  fprintf ( 1, '\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '      H             A           T(H,A)       Exact\n' );
  fprintf ( 1, '\n' );

  n_data = 0;

  while ( 1 )

    [ n_data, h, a, t ] = owen_values ( n_data );

    if ( n_data <= 0 )
      break;
    end

    t2 = tfn ( h, a );

    fprintf ( 1, '  %14f  %14f  %14f  %14f\n', h, a, t2, t );

  end
