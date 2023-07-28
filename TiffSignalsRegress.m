function TiffSignalsRegress(mouse, date, runs, varargin)
%TiffSignalsRegress regresses signal to xyshifts
% TiffSignalsRegress(mouse, date, runs, varargin)


%% Parse
% Whether or not axons was used is determined from icaguidata.pars
p = inputParser;
addOptional(p, 'pmt', 1);  % PMT to use for extraction
addOptional(p, 'server', []);  % Which server to analyss from
addOptional(p, 'optotune', []); % Optotune
addOptional(p, 'force', false);
addOptional(p, 'sigsuffix', '');

% Threshold
addOptional(p, 'input', 'raw'); % Min number of pixels for a cell (after binning)

parse(p, varargin{:});
p = p.Results;

% IO
% Find signals file
if isempty(p.optotune)
    sigpath = sbxPath(mouse, date, runs, 'signals', 'server', p.server);
else
    sigpath = sbxPath(mouse, date, runs, 'OTsig', 'server', p.server, 'pmt', p.pmt, 'optotune', p.optotune);
end
% Multi signal path
if ~isempty(p.sigsuffix)
    sl = length(p.sigsuffix);
    [fp_temp, fn_temp, ext_temp] = fileparts(sigpath);
    
    if ~strcmp(fn_temp(end-sl+1:end), p.sigsuffix)
        fn_alt = sprintf('%s_%s%s', fn_temp, p.sigsuffix, ext_temp);
        sigpath = fullfile(fp_temp, fn_alt);
    end
end

sigstruct = load(sigpath, '-mat');


if isempty(sigstruct.cellsort)
    return;
end

% Check previous reg
if isfield(sigstruct.cellsort(1).timecourse, [p.input, 'regr']) && ~p.force
    fprintf('%s_%s_run%i: Regression already exists, skipping.\n', mouse, date, runs);
    return
end

% Find the shift
if isempty(p.optotune)
    [fp, ~] = fileparts(sigpath);
    fp = fp(1 : find(fp == '\',1,'last'));
    shiftpath = dir(fullfile(fp, 'xyreg', sprintf('%s_%s_%03d*.mat', mouse, date, runs)));
    shiftpath = fullfile(shiftpath.folder, shiftpath.name);
    xyshifts = load(shiftpath, '-mat');
    xyshifts = xyshifts.xy_shifts(:,3:4);
else
    [fp, ~] = fileparts(sigpath);
    fp = fp(1 : find(fp == '\',1,'last'));
    shiftpath = dir(fullfile(fp, 'xyreg', sprintf('%s_%s_%03d_OT%i*.mat', mouse, date, runs, p.optotune)));
    shiftpath = fullfile(shiftpath.folder, shiftpath.name);
    xyshifts = load(shiftpath, '-mat');
    xyshifts = xyshifts.xy_shifts(:,3:4);
end


%% Regress
% new field
newfield = [p.input, 'regr'];
for i = 1 : length(sigstruct.cellsort)
    % data vector
    vec = sigstruct.cellsort(i).timecourse.(p.input);
    
    % weights
    wt = glmfit(xyshifts, vec', 'normal');
    
    % New vec
    vec2 = vec -  wt(2:3)' * xyshifts';
    
    % Save
    sigstruct.cellsort(i).timecourse.(newfield) = vec2;
end

% Save trans
sigstruct.xyshifts = xyshifts;

%% Save
% Save file
save(sigpath, '-struct', 'sigstruct', '-v7.3');
disp('Regression done.')
end