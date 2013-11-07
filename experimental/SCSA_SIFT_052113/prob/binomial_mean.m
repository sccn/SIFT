function mean = binomial_mean ( a, b )

%% BINOMIAL_MEAN returns the mean of the Binomial PDF.
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
%    Input, integer A, the number of trials.
%    1 <= A.
%
%    Input, real B, the probability of success on one trial.
%    0.0 <= B <= 1.0.
%
%    Output, real MEAN, the expected value of the number of
%    successes in A trials.
%
  mean = a * b;
