function cdf = poisson_cdf ( x, a )

%% POISSON_CDF evaluates the Poisson CDF.
%
%  Definition:
%
%    CDF(X,A) is the probability that the number of events observed
%    in a unit time period will be no greater than X, given that the
%    expected number of events in a unit time period is A.
%
%  Modified:
%
%    19 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X, the argument of the CDF.
%    0 <= X.
%
%    Input, real A, the parameter of the PDF.
%    0.0 < A.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x < 0 )

    cdf = 0.0;

  else

    new = exp ( -a );
    sum2 = new;

    for i = 1 : x
      last = new;
      new = last * a / i;
      sum2 = sum2 + new;
    end

    cdf = sum2;

  end
