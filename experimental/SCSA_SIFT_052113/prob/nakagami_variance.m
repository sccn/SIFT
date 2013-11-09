function variance = nakagami_variance ( a, b, c )

%% NAKAGAMI_VARIANCE returns the variance of the Nakagami PDF.
%
%  Modified:
%
%    17 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, C, the parameters of the PDF.
%    0.0 < B
%    0.0 < C
%
%    Output, real VARIANCE, the variance of the PDF.
%
  t1 = gamma ( c + 0.5 );
  t2 = gamma ( c );

  variance = b * b * ( 1.0 - t1 * t1 / ( c * t2 * t2 ) );
