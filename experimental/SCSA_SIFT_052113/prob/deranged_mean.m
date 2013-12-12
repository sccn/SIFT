function mean = deranged_mean ( a )

%% DERANGED_MEAN returns the mean of the Deranged CDF.
%
%  Discussion:
%
%    The mean is computed by straightforward summation.
%
%  Modified:
%
%    07 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer A, the number of items.
%    1 <= A.
%
%    Output, real MEAN, the mean of the PDF.
%
  mean = 0.0;

  for x = 1 : a
    pdf = deranged_pdf ( x, a );
    mean = mean + pdf * x;
  end
