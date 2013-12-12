%% TEST107 tests LORENTZ_MEAN;
%% TEST107 tests LORENTZ_SAMPLE;
%% TEST107 tests LORENTZ_VARIANCE.
%
  clear

  nsample = 1000;
  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST107\n' );
  fprintf ( 1, '  For the Lorentz PDF:\n' );
  fprintf ( 1, '  LORENTZ_MEAN computes the mean;\n' );
  fprintf ( 1, '  LORENTZ_VARIANCE computes the variance;\n' );
  fprintf ( 1, '  LORENTZ_SAMPLE samples.\n' );

  mean = lorentz_mean ( 'DUMMY' );
  variance = lorentz_variance ( 'DUMMY' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF mean =                    %14f\n', mean );
  fprintf ( 1, '  PDF variance =                %14f\n', variance );

  for i = 1 : nsample
    [ x(i), seed ] = lorentz_sample ( seed );
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
