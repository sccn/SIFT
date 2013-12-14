function mean = von_mises_mean ( a, b )

%% VON_MISES_MEAN returns the mean of the von Mises PDF.
%
%  Modified:
%
%    18 September 2004
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
%    Output, real MEAN, the mean of the PDF.
%
  mean = a;
