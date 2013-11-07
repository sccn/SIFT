%% TEST118 tests NORMAL_MEAN;
%% TEST118 tests NORMAL_SAMPLE;
%% TEST118 tests NORMAL_VARIANCE.
%
  clear

  nsample = 1000;
  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST118\n' );
  fprintf ( 1, '  For the Normal PDF:\n' );
  fprintf ( 1, '  NORMAL_MEAN computes the mean;\n' );
  fprintf ( 1, '  NORMAL_SAMPLE samples;\n' );
  fprintf ( 1, '  NORMAL_VARIANCE returns the variance.\n' );

  a = 100.0;
  b = 15.0;

  check = normal_check ( a, b );

  if ( ~check );
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TEST118 - Fatal error!\n' );
    fprintf ( 1, '  The parameters are not legal.\n' );
    return
  end

  mean = normal_mean ( a, b );
  variance = normal_variance ( a, b );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter A =     %14f\n', a );
  fprintf ( 1, '  PDF parameter B =     %14f\n', b );
  fprintf ( 1, '  PDF mean =            %14f\n', mean );
  fprintf ( 1, '  PDF variance =        %14f\n', variance );

  for i = 1 : nsample
    [ x(i), seed ] = normal_sample ( a, b, seed );
  end

  mean = r8vec_mean ( nsample, x );
  variance = r8vec_variance ( nsample, x );
  xmax = max ( x(1:nsample) );
  xmin = min ( x(1:nsample) );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Sample size =     %6d\n', nsample );
  fprintf ( 1, '  Sample mean =     %14f\n', mean );
  fprintf ( 1, '  Sample variance = %14f\n', variance );
  fprintf ( 1, '  Sample maximum =  %14f\n', xmax );
  fprintf ( 1, '  Sample minimum =  %14f\n', xmin );
