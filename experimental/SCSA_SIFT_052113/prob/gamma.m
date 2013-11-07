function value = gamma ( x )

%% GAMMA calculates the Gamma function for a real argument X.
%
%  Definition:
%
%    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
%
%  Recursion:
%
%    GAMMA(X+1) = X * GAMMA(X)
%
%  Special values:
%
%    GAMMA(0.5) = SQRT(PI)
%    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
%
%  Discussion:
%
%    Computation is based on an algorithm outlined in reference 1.
%    The program uses rational functions that approximate the GAMMA
%    function to at least 20 significant decimal digits.  Coefficients
%    for the approximation over the interval (1,2) are unpublished.
%    Those for the approximation for X .GE. 12 are from reference 2.
%    The accuracy achieved depends on the arithmetic system, the
%    compiler, the intrinsic functions, and proper selection of the
%    machine dependent constants.
%
%  Machine dependent constants:
%
%    BETA: radix for the floating-point representation.
%    MAXEXP: the smallest positive power of BETA that overflows.
%    XBIG: the largest argument for which GAMMA(X) is representable
%      in the machine, i.e., the solution to the equation
%      GAMMA(XBIG) = BETA**MAXEXP.
%    XMININ: the smallest positive floating-point number such that
%      1/XMININ is machine representable.
%
%    Approximate values for some important machines are:
%
%                               BETA       MAXEXP        XBIG
%
%    CRAY-1         (S.P.)        2         8191        966.961
%    Cyber 180/855
%      under NOS    (S.P.)        2         1070        177.803
%    IEEE (IBM/XT,
%      SUN, etc.)   (S.P.)        2          128        35.040
%    IEEE (IBM/XT,
%      SUN, etc.)   (D.P.)        2         1024        171.624
%    IBM 3033       (D.P.)       16           63        57.574
%    VAX D-Format   (D.P.)        2          127        34.844
%    VAX G-Format   (D.P.)        2         1023        171.489
%
%                              XMININ
%
%    CRAY-1         (S.P.)   1.84D-2466
%    Cyber 180/855
%      under NOS    (S.P.)   3.14D-294
%    IEEE (IBM/XT,
%      SUN, etc.)   (S.P.)   1.18D-38
%    IEEE (IBM/XT,
%      SUN, etc.)   (D.P.)   2.23D-308
%    IBM 3033       (D.P.)   1.39D-76
%    VAX D-Format   (D.P.)   5.88D-39
%    VAX G-Format   (D.P.)   1.12D-308
%
%  Modified:
%
%    08 October 2004
%
%  Author:
%
%    W. J. Cody and L. Stoltz,
%    Applied Mathematics Division,
%    Argonne National Laboratory,
%    Argonne, Illinois, 60439.
%
%  Reference:
%
%    W J Cody,
%    "An Overview of Software Development for Special Functions",
%    Lecture Notes in Mathematics, 506,
%    Numerical Analysis Dundee, 1975,
%    G. A. Watson (ed.),
%    Springer Verlag, Berlin, 1976.
%
%    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
%    Computer Approximations,
%    Wiley, 1968.
%
%  Parameters:
%
%    Input, real X, the argument of the function.
%
%    Output, real GAMMA, the value of the function.
%    The computation is believed to be free of underflow and overflow.
%
  c = [ ...
    -1.910444077728E-03, ...
     8.4171387781295E-04, ...
    -5.952379913043012E-04, ...
     7.93650793500350248E-04, ...
    -2.777777777777681622553E-03, ...
     8.333333333333333331554247E-02, ...
     5.7083835261E-03 ];

  p = [ ...
    -1.71618513886549492533811E+00, ...
     2.47656508055759199108314E+01, ...
    -3.79804256470945635097577E+02, ...
     6.29331155312818442661052E+02, ...
     8.66966202790413211295064E+02, ...
    -3.14512729688483675254357E+04, ...
    -3.61444134186911729807069E+04, ...
     6.64561438202405440627855E+04 ];

  q = [ ...
    -3.08402300119738975254353E+01, ...
     3.15350626979604161529144E+02, ...
    -1.01515636749021914166146E+03, ...
    -3.10777167157231109440444E+03, ...
     2.25381184209801510330112E+04, ...
     4.75584627752788110767815E+03, ...
    -1.34659959864969306392456E+05, ...
    -1.15132259675553483497211E+05 ];
  sqrtpi = 0.9189385332046727417803297;
  xbig = 35.040;
  xminin = 1.18E-38;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
%
%  Argument is negative.
%
  if ( y <= 0.0 )

    y = -x;
    y1 = floor ( y );
    value = y - y1;

    if ( value ~= 0.0 )

      if ( y1 ~= floor ( y1 * 0.5 ) * 2.0 )
        parity = 1;
      end

      fact = - pi / sin ( pi * value );
      y = y + 1.0;

    else

      value = r8_huge ( value );
      return

    end

  end
%
%  Argument < EPS
%
  if ( y < r8_epsilon ( y ) )

    if ( xminin <= y )
      value = 1.0 / y;
    else
      value = r8_huge ( value );
      return
    end

  elseif ( y < 12.0 )

    y1 = y;
%
%  0.0 < argument < 1.0
%
    if ( y < 1.0 )

      z = y;
      y = y + 1.0;
%
%  1.0 < argument < 12.0, reduce argument if necessary.
%
    else
      n = floor ( y ) - 1;
      y = y - n;
      z = y - 1.0;
    end
%
%  Evaluate approximation for 1.0 < argument < 2.0.
%
    xnum = 0.0;
    xden = 1.0;
    for i = 1 : 8
      xnum = ( xnum + p(i) ) * z;
      xden = xden * z + q(i);
    end

    value = xnum / xden + 1.0;
%
%  Adjust result for case  0.0 < argument < 1.0.
%
    if ( y1 < y )

      value = value / y1;
%
%  Adjust result for case  2.0 < argument < 12.0.
%
    elseif ( y < y1 )

      for i = 1 : n
        value = value * y;
        y = y + 1.0;
      end

    end
%
%  Evaluate for 12 <= argument.
%
  else

    if ( y <= xbig )

      ysq = y * y;
      sum2 = c(7);
      for i = 1 : 6
        sum2 = sum2 / ysq + c(i);
      end
      sum2 = sum2 / y - y + sqrtpi;
      sum2 = sum2 + ( y - 0.5 ) * log ( y );
      value = exp ( sum2 );

    else

      value = r8_huge ( 'DUMMY' );
      return

    end

  end
%
%  Final adjustments and return.
%
  if ( parity )
    value = -value;
  end

  if ( fact ~= 1.0 )
    value = fact / value;
  end
