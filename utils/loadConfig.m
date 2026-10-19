function cfg = loadConfig(file)
%LOADCONFIG Load a configuration struct from a .m function file or a .json file.
%
%   cfg = loadConfig('configs/cfg_fig1.m')    calls the function cfg_fig1()
%   cfg = loadConfig('configs/myrun.json')    decodes the JSON file
%
%   A .m configuration must be a function with no inputs that returns a
%   struct, and its folder must be on the MATLAB path.

[~, name, ext] = fileparts(file);

switch lower(ext)
    case '.m'
        cfg = feval(name);
    case '.json'
        cfg = jsondecode(fileread(file));
    otherwise
        error('loadConfig:unknownType', ...
              'Unknown configuration file type "%s". Use .m or .json.', ext);
end

if ~isstruct(cfg) || ~isscalar(cfg)
    error('loadConfig:notStruct', ...
          'Configuration "%s" did not produce a single struct.', file);
end
end
