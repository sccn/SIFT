%% TEST008 tests BENFORD_PDF.
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST008\n' );
  fprintf ( 1, '  For the Benford PDF:\n' );
  fprintf ( 1, '  BENFORD_PDF evaluates the PDF.\n' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  N    PDF(N)\n' );
  fprintf ( 1, '\n' );

  for n = 1 : 19

    pdf = benford_pdf ( n );
    fprintf ( 1, '  %6d  %14f\n', n, pdf );

  end
