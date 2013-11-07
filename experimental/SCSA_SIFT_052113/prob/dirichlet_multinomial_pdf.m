function pdf = dirichlet_multinomial_pdf ( x, a, b, c )

%% DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
%
%  Formula:
%
%    PDF(X)(A,B,C) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
%      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
%
%    where:
%
%      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
%      C_Sum = Sum ( 1 <= I <= B ) C(I)
%
%  Modified:
%
%    19 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Kenneth Lange,
%    Mathematical and Statistical Methods for Genetic Analysis,
%    Springer, 1997, page 45.
%
%  Parameters:
%
%    Input, integer X(B); X(I) counts the number of occurrences of
%    outcome I, out of the total of A trials.
%
%    Input, integer A, the total number of trials.
%
%    Input, integer B, the number of different possible outcomes on
%    one trial.
%
%    Input, real C(B); C(I) is the Dirichlet parameter associated
%    with outcome I.
%
%    Output, real PDF, the value of the Dirichlet multinomial PDF.
%
  c_sum = sum ( c(1:b) );

  pdf_log = -gamma_log ( c_sum + a ) + gamma_log ( c_sum ) + gamma_log ( a + 1 );

  for i = 1 : b
    pdf_log = pdf_log + gamma_log ( c(i) + x(i) ) - gamma_log ( c(i) ) ...
      - gamma_log ( x(i) + 1 );
  end

  pdf = exp ( pdf_log );
