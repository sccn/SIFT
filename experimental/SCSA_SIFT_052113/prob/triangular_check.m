function check = triangular_check ( a, b )

%% TRIANGULAR_CHECK checks the parameters of the Triangular CDF.
%
%  Modified:
%
%    22 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    A < B.
%
%    Output, logical CHECK, is true if the parameters are legal.
%
  if ( b <= a )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TRIANGULAR_CHECK - Fatal error!\n' );
    fprintf ( 1, '  B <= A.\n' );
    check = 0;
    return
  end

  check = 1;

