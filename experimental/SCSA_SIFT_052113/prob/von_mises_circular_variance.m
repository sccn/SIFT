function value = von_mises_circular_variance ( a, b )

%% VON_MISES_CIRCULAR_VARIANCE returns the circular variance of the von Mises PDF.
%
%  Modified:
%
%    02 December 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, B, the parameters of the PDF.
%    -PI <= A <= PI,
%    0.0 < B.
%
%    Output, real VALUE, the circular variance of the PDF.
%
  value = 1.0 - bessel_i1 ( b ) / bessel_i0 ( b );
