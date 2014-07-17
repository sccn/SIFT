function cdf = maxwell_cdf ( x, a )

%% MAXWELL_CDF evaluates the Maxwell CDF.
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
%    Input, real X, the argument of the PDF.
%    0.0 <= X
%
%    Input, real A, the parameter of the PDF.
%    0 < A.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= 0.0 )

    cdf = 0.0;

  else

    x2 = x / a;
    p2 = 1.5;

    cdf = gamma_inc ( p2, x2 );

  end
