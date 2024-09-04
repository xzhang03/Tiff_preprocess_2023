function Tiffxyreg(mouse, date, run, varargin)
%Tiffxyreg applies xyreg to tiff movies
% TiffXYReg(mouse, date, run, varargin)

%% Parse inputs
p = inputParser;

addOptional(p, 'server', []);  % Which server to call, if not the current server
addOptional(p, 'pmt', []);  % Which PMT, in the case of sbxreg and sbxclean
addOptional(p, 'input', 'OTtiff'); % Input
addOptional(p, 'optotune', []); % Which optotune, used for optotune-split tiffs.
addOptional(p, 'filepath', ''); % Directly feed in filepath

% Focus area (area used to calculate shifts)
addOptional(p, 'autofocus', true); % If true, there is no manual selection of focus area
addOptional(p, 'autofocuspos', [0.25 0.75 0.25 0.75]);% On a 0 to 1 scale as [Row_start, Row_end, Column_start, Column_end]

% Edges and downsample
addOptional(p, 'binxy', [], @isnumeric);  % The maximum cross-correlation to allow between overlapping ROIs, combined with overlap, default 2
addOptional(p, 'edges', []);  % The edges of the image to be removed before ROI extraction. Will be set to sbxRemoveEdges if empty

% Ref
addOptional(p, 'reuseref', false); % Reuse reference if can find it (must be the same size).
addOptional(p, 'refsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
addOptional(p, 'refoffset', 100, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames

% Gaussian sizes
addOptional(p, 'gausssize', [8 30]);

% Parallel processing
addOptional(p, 'useparfor', true);
addOptional(p, 'chunksize', 1000); % Chunk size for parallel processing. Decrease if RAM is an issue

% Iterative registration
addOptional(p, 'nitr', 1);

if length(varargin) == 1 && iscell(varargin{1}), varargin = varargin{1}; end
parse(p, varargin{:});
p = p.Results;

%% IO
if ~isempty(p.filepath)
    datapath = p.filepath;
elseif isempty(p.optotune)
    % No optotune
    datapath = sbxPath(mouse,date,run,'sbx','pmt', p.pmt, 'server', p.server);
    if ~strcmpi(p.input, 'sbx')
        [fp, fn, ~] = fileparts(datapath);
        datapath = fullfile(fp, [fn, '.tif']);

        if ~exist(datapath, 'file')
            error('No tiff was found.')
        end
    end
else
    % Optotune split
    datapath = sbxPath(mouse, date, run, p.input, 'pmt', p.pmt, 'optotune', p.optotune, 'server', p.server);
    if ~exist(datapath, 'file')
        error('No tiff was found.')
    end
end

% Read
switch p.input
    case 'sbx'
        tic;
        fprintf('Reading sbx... ');
        im = sbxReadPMT(datapath, 0, -1, 0, p.optotune);
        fprintf('Done. Elapsed time (s) = %i\n', round(toc));
    otherwise
        im = readtiff(datapath);
end
fprintf('Initializing...\n')
tic;
type = class(im);

% Output filnames
[fp, fn] = fileparts(datapath);
fn_out = fullfile(fp, sprintf('%s_xyreg-%i.tif', fn, p.pmt));
ref_out = fullfile(fileparts(fp), 'xyreg', sprintf('%s_xyregref-%i.tif', fn, p.pmt));
fn_shifts = fullfile(fileparts(fp), 'xyreg', sprintf('%s_xyreg-%i.mat', fn, p.pmt));

% Make xyreg folder if necessary
if ~exist(fullfile(fileparts(fp), 'xyreg'), 'dir')
    mkdir(fullfile(fileparts(fp), 'xyreg'));
end

hwait = waitbar(0, 'Processing.');
for itr = 1 : p.nitr
    %% Preprocess and ref
    % Apply edges
    im2 = single(im(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :));
    
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Binning.', itr, p.nitr));
    % downsample
    if p.binxy > 1
        im2bin = binxy(im2);
    else
        im2bin = im2;
    end

    % ref
    if exist(ref_out, 'file') && p.reuseref
        fprintf('Reusing reference...\n');
        waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Reusing reference', itr, p.nitr));
        imref = readtiff(ref_out);
    else
        imref = mean(im2bin(:,:,p.refoffset : p.refoffset + p.refsize - 1), 3);
        fprintf('Making reference...\n');
        waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Making reference', itr, p.nitr));
        writetiff(imref, ref_out);
    end

    % Get a focus area for calculating shifts
    if ~p.autofocus
        imshow(imref,[]);
        h = imrect;
        regfocus = wait(h);
        regfocus = round(regfocus);
        close(gcf)

        % Get coordinates
        regfocus2 = [regfocus(2), regfocus(2) + regfocus(4) - 1, regfocus(1), regfocus(1) + regfocus(3) - 1];
    else
        % Get focus
        sz = size(imref);
        regfocus2 = round([sz(1), sz(1), sz(2), sz(2)] .* p.autofocuspos);
    end

    % Crop reference
    imref = imref(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4));

    % Normalize reference
    imref = medfilt2(imref, [2 2], 'symmetric');
    f_prime = imref - imgaussfilt(imref, p.gausssize(1));
    imref = f_prime ./ (imgaussfilt(f_prime.^2, p.gausssize(2)).^(1/2));
    imref(isnan(imref)) = 0;

    t = round(toc);
    fprintf('Initilization done. Elapsed time = %i seconds.\n', t);
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Initilization done. Elapsed time = %i seconds.', itr, p.nitr, t));



    %% Parallel register
    % Number of frames
    nframes = size(im, 3);

    % Parpool
    if p.useparfor
        if isempty(gcp('nocreate'))
            parpool();
        end
    end

    % Parallel processing
    nchunks = ceil(nframes/p.chunksize);

    % Display
    tnow = datestr(datetime('now'), 13);
    fprintf('Parallel registration starting at %s.\n', tnow)
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Parallel registration starting at %s.', itr, p.nitr, tnow));

    
    % Initialize shifts
    xy_shifts = cell(nchunks, 1);
    tic;
    parfor i = 1 : nchunks
        % Call core function
        xy_shifts{i} = Tiffxyreg_core(im2bin(regfocus2(1) : regfocus2(2), regfocus2(3) : regfocus2(4), (i-1) * p.chunksize + 1:...
            min(i * p.chunksize, nframes)), imref, p.gausssize);
    end


    xy_shifts = cell2mat(xy_shifts);
    t = round(toc);    
    fprintf('\nParallel registration done. Elapsed time = %i seconds.\n', t);
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Parallel registration done. Elapsed time = %i seconds.', itr, p.nitr, t));
    
    % Modify the shifts according to binning
    if p.binxy > 1
        xy_shifts(:,3:4) = xy_shifts(:,3:4) * p.binxy;
    end

    %% Apply shifts
    % Apply xyshifts
    tnow = datestr(datetime('now'), 13);
    fprintf('Registration application starting at %s.\n', tnow)
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Registration application starting at %s.', itr, p.nitr, tnow));
    
    tic;
    [~, imreg] = stackRegisterMA_RR(im2, [], [], xy_shifts);

    % Apply edges
    im(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :) = cast(imreg, type);

    % Apply class
    t = round(toc);    
    fprintf('\nRegistration applied. Elapsed time = %i seconds.\n', t);
    waitbar(itr/p.nitr, hwait, sprintf('%i/%i: Registration applied. Elapsed time = %i seconds.', itr, p.nitr, t));

end
close(hwait)

%% Save
% Shifts
save(fn_shifts, '-v7.3', 'xy_shifts', 'p', 'imref', 'nframes');
fprintf('Shifts saved.\n')

% Image
writetiff(im, fn_out);
fprintf('Done.\n')
end