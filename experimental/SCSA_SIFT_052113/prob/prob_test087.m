%% TEST087 tests GUMBEL_MEAN;
%% TEST087 tests GUMBEL_SAMPLE;
%% TEST087 tests GUMBEL_VARIANCE.
%
  clear

  nsample = 1000;
  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST087\n' );
  fprintf ( 1, '  For the Gumbel PDF:\n' );
  fprintf ( 1, '  GUMBEL_MEAN computes the mean;\n' );
  fprintf ( 1, '  GUMBEL_SAMPLE samples;\n' );
  fprintf ( 1, '  GUMBEL_VARIANCE computes the variance.\n' );

  mean = gumbel_mean ( 'DUMMY' );

  variance = gumbel_variance ( 'DUMMY' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF mean      =               %14f\n', mean );
  fprintf ( 1, '  PDF variance =                %14f\n', variance );

  for i = 1 : nsample
    [ x(i), seed ] = gumbel_sample ( seed );
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
