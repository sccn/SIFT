function cdf = triangle_cdf ( x, a, b, c )

%% TRIANGLE_CDF evaluates the Triangle CDF.
%
%  Modified:
%
%    21 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the CDF.
%
%    Input, real A, B, C, the parameters of the PDF.
%    A <= B <= C and A < C.
%
%    Output, real CDF, the value of the CDF.
%
  if ( x <= a )

    cdf = 0.0;

  elseif ( x <= b )

    if ( a == b )
      cdf = 0.0;
    else
      cdf = ( x - a ) * ( x - a ) / ( b - a ) / ( c - a );
    end

  elseif ( x <= c )

    cdf = ( b - a ) / ( c - a ) + ( 2.0 * c - b - x ) * ( x - b ) / ( c - b ) / ( c - a );

  else

    cdf = 1.0;

  end

