


if isfield(opts,'FileNameLoad')
    opts.WaveOneFileName = opts.FileNameLoad;
    opts=rmfield(opts,'FileNameLoad')
end

if isfield(opts,'FileNameLoad2')
    opts.WaveTwoFileName = opts.FileNameLoad2;
    opts=rmfield(opts,'FileNameLoad2')
end

if isfield(opts,'FileNameIndex')
       opts.WaveOneFileIndex = opts.FileNameIndex;
       opts=rmfield(opts,'FileNameIndex')
end

if isfield(opts,'FileNameIndex2')
       opts.WaveTwoFileIndex = opts.FileNameIndex2;
       opts=rmfield(opts,'FileNameIndex2')
end



