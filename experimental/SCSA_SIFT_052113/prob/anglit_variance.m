function variance = anglit_variance ( dummy )

%% ANGLIT_VARIANCE returns the variance of the Anglit PDF.
%
%  Discussion:
%
%    Variance =
%      Integral ( -PI/4 <= X <= PI/4 ) X**2 * SIN ( 2 * X + PI / 2 )
%
%    Antiderivative =
%      0.5 * X * SIN ( 2 * X + PI / 2 )
%      + ( 0.25 - 0.5 * X**2 ) * COS ( 2 * X + PI / 2 )
%
%  Modified:
%
%    01 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = 0.0625 * pi^2 - 0.5;
