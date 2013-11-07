%% TEST002 tests ANGLE_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST002\n' );
  fprintf ( 1, '  For the ANGLE PDF:\n' );
  fprintf ( 1, '  ANGLE_PDF evaluates the PDF;\n' );

  n = 5;
  x = 0.50;

  pdf = angle_pdf ( x, n );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  PDF parameter N = %6d\n', n );
  fprintf ( 1, '  PDF argument X =  %14f\n', x );
  fprintf ( 1, '  PDF value =       %14f\n', pdf );
