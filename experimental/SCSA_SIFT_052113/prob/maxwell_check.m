function check = maxwell_check ( a )

%% MAXWELL_CHECK checks the parameters of the Maxwell CDF.
%
%  Modified:
%
%    16 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, the parameter of the PDF.
%    0 < A.
%
%    Output, logical CHECK, is true if the parameters are legal.
%
  if ( a <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'MAXWELL_CHECK - Fatal error!\n' );
    fprintf ( 1, '  A <= 0.0.\n' );
    check = 0;
    return
  end

  check = 1;
