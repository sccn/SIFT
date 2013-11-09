function variance = i4vec_variance ( n, x )

%% I4VEC_VARIANCE returns the variance of an I4VEC.
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
%    Input, integer X(N), the vector whose variance is desired.
%
%    Output, real VARIANCE, the variance of the vector entries.
%
  mean = i4vec_mean ( n, x );

  variance = sum ( ( x(1:n) - mean ).^2 );

  if ( 1 < n )
    variance = variance / ( n - 1 );
  else
    variance = 0.0;
  end 

