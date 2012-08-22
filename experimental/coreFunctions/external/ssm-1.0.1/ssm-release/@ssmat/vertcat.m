function A = vertcat(varargin)

%@SSMAT/VERTCAT Vertical concatenation of SSMAT objects.

% (c) 2006-2007 Jyh-Ying Peng �^����
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

n       = 1;
mat     = cell(1, nargin);
mmask   = cell(1, nargin);
dmmask  = cell(1, nargin);
dvec    = cell(1, nargin);
dvmask  = cell(1, nargin);
for i = 1 : nargin
    if isa(varargin{i}, 'ssmat'), m = varargin{i};
    elseif isnumeric(varargin{i}) && ndims(varargin{i}) <= 3, m = ssmat(varargin{i});
    else error('ssm:ssmat:vertcat:UnableToConvert', ['Input ' int2str(i) ' cannot be converted to SSMAT class.']);
    end
    if n == 1, n = size(m.dvec, 2);
    elseif size(m.dvec, 2) ~= 1 && n ~= size(m.dvec, 2), error('ssm:ssmat:vertcat:InputError', 'vertcat is not defined for SSMAT with different time durations.');
    end
    if size(m, 2) ~= 1, error('ssm:ssmat:vertcat:InputError', 'vertcat is only defined for column vector SSMATs.'); end
    mat{i}      = m.mat;
    if isempty(m.mmask), mmask{i} = false(size(m.mat)); else mmask{i} = m.mmask; end
    if isempty(m.dmmask), dmmask{i} = false(size(m.mat)); else dmmask{i} = m.dmmask; end
    dvec{i}     = m.dvec;
    if isempty(m.dvmask), dvmask{i} = false(size(m.dvec, 1), 1); else dvmask{i} = m.dvmask; end
end
A.mat      = vertcat(mat{:});
A.mmask    = vertcat(mmask{:});
A.dmmask   = vertcat(dmmask{:});
A.dvec     = vertcat(dvec{:});
A.dvmask   = vertcat(dvmask{:});

%% Eliminate mask degeneracy %%
if ~isempty(A.mmask) && all(~A.mmask(:)), A.mmask = []; end
if ~isempty(A.dmmask) && all(~A.dmmask(:)), A.dmmask = []; end
if ~isempty(A.dvmask) && all(~A.dvmask(:)), A.dvmask = []; end

%% Register object instance %%
A   = class(A, 'ssmat');

