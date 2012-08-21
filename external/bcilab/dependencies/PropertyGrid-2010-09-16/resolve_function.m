function f = resolve_function(expr)
% resolve a scoped function, which is normally not possible in MATLAB

% thus, we have to use a crazy trick: arg_report gives us access to nested functions
ofs = strfind(expr,'::@');
f = arg_report('handle',str2func(expr(1:ofs-1)),{expr(ofs+3:end)});
f = f{1};
