function value = r8poly_value ( n, a, x )

%% R8POLY_VALUE evaluates a double precision polynomial.
%
%  Discussion:
%
%    For sanity's sake, the value of N indicates the NUMBER of 
%    coefficients, or more precisely, the ORDER of the polynomial,
%    rather than the DEGREE of the polynomial.  The two quantities
%    differ by 1, but cause a great deal of confusion.
%
%    Given N and A, the form of the polynomial is:
%
%      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
%
%  Modified:
%
%    12 August 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the order of the polynomial.
%
%    Input, real A(1:N), the coefficients of the polynomial.
%    A(1) is the constant term.
%
%    Input, real X, the point at which the polynomial is to be evaluated.
%
%    Output, real VALUE, the value of the polynomial at X.
%
  value = 0.0;
  for i = n : -1 : 1
    value = value * x + a(i);
  end
