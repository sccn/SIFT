function mean = bernoulli_mean ( a )

%% BERNOULLI_MEAN returns the mean of the Bernoulli PDF.
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
%    Input, real A, the probability of success.
%    0.0D+00 <= A <= 1.0.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = a;
