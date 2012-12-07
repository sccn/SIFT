function res = unassigned
% special keyword (constant function) that can be used in arg*() specifications, 
% in particular for the default values of arguments that do not show up in the function's workspace (or the GUI) unless explicitly assigned
res = '__arg_unassigned__';