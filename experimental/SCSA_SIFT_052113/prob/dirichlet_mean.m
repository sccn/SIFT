function mean = dirichlet_mean ( n, a )

%% DIRICHLET_MEAN returns the means of the Dirichlet PDF.
%
%  Modified:
%
%    30 October 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of components.
%
%    Input, real A(N), the probabilities for each component.
%    Each A(I) should be positive.
%
%    Output, real MEAN(N), the means of the PDF.
%
  mean(1:n) = r8vec_unit_sum ( n, a );
