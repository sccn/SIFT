function mean = i4vec_mean ( n, x )

%% I4VEC_MEAN returns the mean of an I4VEC.
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
%    Input, integer N, the number of entries in the vector.
%
%    Input, integer X(N), the vector whose mean is desired.
%
%    Output, real MEAN, the mean, or average,
%    of the vector entries.
%
  mean = sum ( x(1:n) ) / n;
