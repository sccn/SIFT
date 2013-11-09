function variance = bernoulli_variance ( a )

%% BERNOULLI_VARIANCE returns the variance of the Bernoulli PDF.
%
%  Modified:
%
%    02 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real A, the probability of success on one trial.
%    0.0D+00 <= A <= 1.0.
%
%    Output, real VARIANCE, the variance of the PDF.
%
  variance = a * ( 1.0 - a );
