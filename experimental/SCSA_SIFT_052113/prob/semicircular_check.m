function check = semicircular_check ( a, b )

%% SEMICIRCULAR_CHECK checks the parameters of the Semicircular CDF.
%
%  Modified:
%
%    20 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameter of the PDF.
%    0.0 < B.
%
  if ( b <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SEMICIRCULAR_CHECK - Fatal error!\n' );
    fprintf ( 1, '  B <= 0.0\n' );
    check = 0;
    return
  end

  check = 1;
