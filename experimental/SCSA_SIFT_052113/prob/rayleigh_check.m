function check = rayleigh_check ( a )

%% RAYLEIGH_CHECK checks the parameter of the Rayleigh PDF.
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
%    Input, real A, the parameter of the PDF.
%    0.0 < A.
%
%    Output, logical CHECK, is true if the parameters are legal.
%
  if ( a <= 0.0 )
    fprintf ( 1, '\n');
    fprintf ( 1, 'RAYLEIGH_CHECK - Fatal error!\n' );
    fprintf ( 1, '  A <= 0.\n' );
    check = 0;
    return
  end

  check = 1;

