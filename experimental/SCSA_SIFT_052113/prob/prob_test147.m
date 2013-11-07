%% TEST147 tests UNIFORM_01_ORDER_SAMPLE;
%
  clear

  n = 10;
  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST147\n' );
  fprintf ( 1, '  For the Uniform 01 Order PDF:\n' );
  fprintf ( 1, '  UNIFORM_ORDER_SAMPLE samples.\n' );
  fprintf ( 1, '\n' );

  [ x, seed ] = uniform_01_order_sample ( n, seed );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Ordered sample:\n' );
  fprintf ( 1, '\n' );

  for i = 1 : n
    fprintf ( 1, '  %6d  %14f\n', i, x(i) );
  end
