%% TEST011 tests BETA;
%% TEST011 tests GAMMA.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST011\n' );
  fprintf ( 1, '  BETA evaluates the Beta function;\n' );
  fprintf ( 1, '  GAMMA evaluates the Gamma function.\n' );

  a = 2.2;
  b = 3.7;

  beta1 = beta ( a, b );
  beta2 = gamma ( a ) * gamma ( b ) / gamma ( a + b );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Argument A =                   %14f\n', a );
  fprintf ( 1, '  Argument B =                   %14f\n', b );
  fprintf ( 1, '  Beta(A,B) =                    %14f\n', beta1 );
  fprintf ( 1, '  (Expected value = 0.0454 )\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Gamma(A)*Gamma(B)/Gamma(A+B) = %14f\n', beta2 );
