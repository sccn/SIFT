%% TEST005 tests ANGLIT_MEAN;
%% TEST005 tests ANGLIT_SAMPLE;
%% TEST005 tests ANGLIT_VARIANCE.
%
  clear

  nsample = 1000;
  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST005\n' );
  fprintf ( 1, '  For the Anglit PDF:\n' );
  fprintf ( 1, '  ANGLIT_MEAN computes the mean;\n' );
  fprintf ( 1, '  ANGLIT_SAMPLE samples;\n' );
  fprintf ( 1, '  ANGLIT_VARIANCE computes the variance.\n' );

  mean = anglit_mean ( 'DUMMY' );
  variance = anglit_variance ( 'DUMMY' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF mean =     %14f\n', mean );
  fprintf ( 1, '  PDF variance = %14f\n', variance );

  for i = 1 : nsample
    [ x(i), seed ] = anglit_sample ( seed );
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

