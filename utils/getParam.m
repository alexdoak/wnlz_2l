function val = getParam(cfg, name, promptFcn, isValid)
%GETPARAM Take a parameter from a configuration struct, or prompt for it.
%
%   val = getParam(cfg, name, promptFcn)
%   returns the field NAME of the struct CFG if it exists; otherwise calls
%   promptFcn() to ask the user. NAME may be nested with dots, e.g.
%   'Bathymetry.stepHeight'.
%
%   val = getParam(cfg, name, promptFcn, isValid)
%   also checks a value taken from CFG with the function handle isValid.
%   If the check fails, a warning is printed and the user is prompted.
%
%   Example:
%     H = getParam(cfg, 'H', @() chooseFloatValue("Total depth H? ",0,[],1), ...
%                  @(x) isnumeric(x) && isscalar(x) && x > 0);

if nargin < 4
    isValid = [];
end

[found, val] = getNestedField(cfg, name);

if found
    if isempty(isValid) || isValid(val)
        fprintf("  %s = %s  (from configuration)\n", name, describeValue(val));
        return
    end
    fprintf("  WARNING: configuration value %s = %s is not valid here. " + ...
            "Please choose it below.\n", name, describeValue(val));
end

val = promptFcn();
end

% ----------------------------------------------------------------------- %
function [found, val] = getNestedField(s, name)
% Follow a dotted field name through nested structs.
found = false;
val = [];
parts = strsplit(name, '.');
for k = 1:numel(parts)
    if ~isstruct(s) || ~isscalar(s) || ~isfield(s, parts{k})
        return
    end
    s = s.(parts{k});
end
found = true;
val = s;
end

% ----------------------------------------------------------------------- %
function str = describeValue(val)
% Short printable description of a value, for the "(from configuration)" line.
if isnumeric(val) || islogical(val)
    str = mat2str(val);
elseif ischar(val) || isstring(val)
    str = sprintf('"%s"', val);
else
    str = ['<' class(val) '>'];
end
end
