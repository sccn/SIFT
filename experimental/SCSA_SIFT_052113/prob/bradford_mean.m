function mean = bradford_mean ( a, b, c )

%% BRADFORD_MEAN returns the mean of the Bradford PDF.
%
%  Modified:
%
%    03 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, C, the parameters of the PDF.
%    A < B,
%    0.0 < C.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = ( c * ( b - a ) + log ( c + 1.0 ) * ( a * ( c + 1.0 ) - b ) ) ...
    / ( c * log ( c + 1.0 ) );

