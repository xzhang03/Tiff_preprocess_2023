function TiffSignalsAddSomaROIs(mouse, date, runs, varargin)
%TiffSignalsCellpose makes signal file from cellpose segmentations. This is
%a tiff-based code.
%   run simplifycellsort

% Parse
% Whether or not axons was used is determined from icaguidata.pars
p = inputParser;
addOptional(p, 'pmt', 1);  % PMT to use for extraction
addOptional(p, 'server', []);  % Which server to analyss from
addOptional(p, 'optotune', []); % Optotune
addOptional(p, 'force', false);  % Overwrite files if they exist
addOptional(p, 'sigsuffix', '');

parse(p, varargin{:});
p = p.Results;

%% IO
% Signals path
if isempty(p.optotune)
    sigpath = sbxPath(mouse, date, runs, 'signals', 'server', 'nasquatch');
    [fp, ~, ~] = fileparts(sigpath);
    sigpath = fullfile(fp,sprintf('%s_%s_%03d.signals', mouse, date, runs));
else
    sigpath = sbxPath(mouse, date, runs, 'OTsig', 'server', p.server, 'optotune', p.optotune);
end

% Multi signal path
if ~isempty(p.sigsuffix)
    [fp_temp, fn_temp, ext_temp] = fileparts(sigpath);
    fn_alt = sprintf('%s_%s%s', fn_temp, p.sigsuffix, ext_temp);
    sigpath = fullfile(fp_temp, fn_alt);
end

% Load
sigstruct = load(sigpath, '-mat');
if isfield(sigstruct, 'movROIs') && ~p.force
    redo = input('ROIs summary already exists, redo? (1 = yes, 0 - no): ');
    if redo ~= 1
        return;
    end
end

%% Calculate
l = length(sigstruct.cellsort);
for i = 1 : l
    switch i
        case 1
            movROIs = sigstruct.cellsort(i).mask > 0;
        otherwise
            movROIs = movROIs + sigstruct.cellsort(i).mask > 0;
    end
end
movROIs = movROIs > 0;

%% Save
save(sigpath, 'movROIs', '-append');

end

