function grinface = TiffSignalsCellpose(mouse, date, runs, varargin)
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
addOptional(p, 'filepath', ''); % Feed tif path directly

% Threshold
addOptional(p, 'minarea', 20); % Min number of pixels for a cell (after binning)
addOptional(p, 'mintm', 1500); % Min tm value to be considered (filter out gross fitting errors)
addOptional(p, 'allowsplit', false); % Allow masks that are split (axons and dendrites)

% IO
addOptional(p, 'movtype', 'OTtiff_demonsreg');  % input type, can be xyreg, sbxreg, or sbx.
addOptional(p, 'binxy', 1);  % Downsample factor used in cellpose
addOptional(p, 'usetrim', true); % Correct for small size differences between movie and segmentation

% Masks
addOptional(p, 'presegsuffix', 'toseg.tif');
addOptional(p, 'masksuffix', '_cp_masks.png');

% Multi signal files (e.g., soma vs dendrite). 
addOptional(p, 'sigsuffix', ''); % Leave empty if not
addOptional(p, 'nosoma', false); % Remove soma overlap

% Neuropil
addOptional(p, 'npsize', [14, 6]); % Neuropil size, avoidance size

% Dff
addOptional(p, 'dffwindow', 32); % Window for dff (in seconds)
addOptional(p, 'freq', []); % fps, needed for dff calculations
addOptional(p, 'percentile', 10); % Percentile for dff calculations

% Edge
addOptional(p, 'edges', []); % Check edges
addOptional(p, 'GRIN', false); % Remove anything that is not within the GRIN surface
addOptional(p, 'reusegrinface', true); % Try to reuse grin face from a previous run
addOptional(p, 'grinface', []); % Pass the face of the GRIN lens here if you don't want manual selection


parse(p, varargin{:});
p = p.Results;

if p.binxy > 1 && ~isempty(p.edges)
    p.edges = round(p.edges/p.binxy);
end

% Alternative signal filename
sigalt = ~isempty(p.sigsuffix);

%% IO
% Tiff path
if ~isempty(p.filepath)
    tiffpath = p.filepath;
else
    tiffpath = sbxPath(mouse, date, runs, p.movtype, 'server', p.server, 'pmt', p.pmt, 'optotune', p.optotune);
end

% Mask path
[fp, ~] = fileparts(tiffpath);

if ~isempty(p.filepath)
    maskpath = dir(sprintf('%s*%s', tiffpath(1:end-4), p.masksuffix));
elseif isempty(p.optotune)
    maskpath = dir(fullfile(fp,[sprintf('%s_%s_%03d*', mouse, date, runs),p.masksuffix]));
else
    maskpath = dir(fullfile(fp,[sprintf('%s_%s_%03d_OT%i*', mouse, date, runs, p.optotune),p.masksuffix]));
end
maskpath = fullfile(maskpath.folder, maskpath.name);

% Presegmentation image path
% meantifpath
if ~isempty(p.filepath)
    meanpath = sprintf('%s_%s', tiffpath(1:end-4), p.presegsuffix);
elseif isempty(p.optotune)
    meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runs));
else
    meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runs,p.optotune));
end
movmean = readtiff(meanpath);

% Signals path
if ~isempty(p.filepath)
    sigpath = sprintf('%s%s', tiffpath(1:end-4), '.mat');
elseif isempty(p.optotune)
    sigpath = fullfile(fp,sprintf('%s_%s_%03d.signals', mouse, date, runs));
else
    sigpath = fullfile(fp,sprintf('%s_%s_%03d_OT%i.signals', mouse, date, runs, p.optotune));
end

% Multi signal path
if sigalt
    [fp_temp, fn_temp, ext_temp] = fileparts(sigpath);
    fn_alt = sprintf('%s_%s%s', fn_temp, p.sigsuffix, ext_temp);
    sigpath_alt = fullfile(fp_temp, fn_alt);
end

% Check exist
if ~sigalt
    % Not alternative signal file
    if exist(sigpath, 'file') && ~p.force
        redo = input('Signal file already exist, redo? (1 = yes, 0 - no): ');
        if redo ~= 1
            grinface = p.grinface;
            return;
        end
    end
else
    % Alternative signal file
    if exist(sigpath_alt, 'file') && ~p.force
        redo = input('Signal file already exist, redo? (1 = yes, 0 - no): ');
        if redo ~= 1
            grinface = p.grinface;
            return;
        end
    end
end

% Get the face of GRIN
if p.GRIN
    if isempty(p.grinface)
        if p.reusegrinface && sigalt && exist(sigpath_alt, 'file')
            loaded = load(sigpath_alt, '-mat', 'psig');
            p.grinface = loaded.psig.grinface;
        elseif p.reusegrinface && exist(sigpath, 'file')
            loaded = load(sigpath, '-mat', 'psig');
            p.grinface = loaded.psig.grinface;
        else
            p.grinface = getpoly(movmean, 'Select the face of the GRIN lens.');
        end
    end
end

% Load ROI summary
if p.nosoma
    somaROIs = load(sigpath, '-mat', 'movROIs');
    somaROIs = somaROIs.movROIs;
    if p.binxy > 1
        somaROIs = binxy(somaROIs);
    end
    somaROIs = somaROIs > 0;
end

% Check dimensions
sampletiff = imread(tiffpath, 1);
try
    masks = readtiff(maskpath);
catch
    masks = imread(maskpath);
end

% Number of masks
nmasks = max(masks(:));

% Continue to check dimensions
if any(size(sampletiff) ~= size(masks) * p.binxy)
    % Check if tolerance is zero
    if ~p.usetrim
        disp('Dimensions do not match')
        grinface = p.grinface;
        return
    end
    
    trimmed = true;
    movtrim = [0 0];
    masktrim = [0 0];
    
    % Check rows
    rowdiff = size(sampletiff,1) - size(masks,1) * p.binxy;
    if abs(rowdiff) >= p.binxy
        disp('The number of rows differ too much.');
    elseif rowdiff > 0
        movtrim(1) = rowdiff;
    elseif rowdiff < 0
        masktrim(1) = -rowdiff;
    end
    
    % Check columns
    coldiff = size(sampletiff,2) - size(masks,2) * p.binxy;
    if abs(coldiff) >= p.binxy
        disp('The number of columns differ too much.');
    elseif coldiff > 0
        movtrim(2) = coldiff;
    elseif coldiff < 0
        masktrim(2) = -coldiff;
    end
else
    trimmed = false;
end

% Output mask as RGB
size(masks)
[fp, fn, ~] = fileparts(maskpath);
if ~sigalt
    % Not alt signal file
    maskrgbpath = fullfile(fp, [fn, '_RGB.png']);
else
    % Alternative signal file
    maskrgbpath = fullfile(fp, [fn, '_', p.sigsuffix, '_RGB.png']);
end
if ~exist(maskrgbpath, 'file')
    % Make rgb
    maskrgb = repmat(imresize(mat2gray(movmean),p.binxy), [1 1 3]);
    ROIs_bin = imresize(mat2gray(masks >= 1),p.binxy);
    maskrgb(:,:,1) = maskrgb(:,:,1) + ROIs_bin;
    
    % Show rgb
    figure
    imshow(maskrgb);
    
    % Label
    for i = 1 : nmasks
        if sum(sum(masks == i)) > 0
            cen = regionprops(masks == i, 'Centroid');
            cen = cen.Centroid;
            cen = cen * p.binxy;
            text(cen(1), cen(2), num2str(i), 'Color', [0 0 0], 'FontSize', 14);
        end
    end
    
    saveas(gcf, maskrgbpath);
    close(gcf)
end

% Read the whole movie
mov = readtiff(tiffpath);

% shrink video (the size here has to match)
fprintf('Binning...');
tic
if p.binxy > 1
    mov = binxy(mov, p.binxy);
end

% Size vec
sz = size(mov);

% In case it doesn't match in size (unlikely)
if any(sz(1:2) ~= size(masks))
    disp('Something went wrong in binning.');
    return
end

% Display
t = round(toc);
fprintf(' Done. Elapsed time = %i seconds.\n', t);

%% Initialize cellsort structure
% Initialize
fprintf('Initializing structure...');
tic
cellsort = struct('mask', [], 'neuropil', [], 'group_number', 1, 'timecourse', [], 'area', [],...
    'centroid', [], 'bbox', []);
cellsort = repmat(cellsort, [max(masks(:)), 1]);

% Dilate mask
masksdilated = imresize(masks > 0, p.binxy);

% Make sure that the mask dimension matches the original movie
if trimmed
    if any(masktrim ~= 0)
        masksdilated = masksdilated(1:end-masktrim(1), 1:end-masktrim(2));
    end

    if movtrim(1) > 0
        masksdilated(end+1 : end+movtrim(1), :) = 0; 
    end

    if movtrim(2) > 0
        masksdilated(:, end+1 : end+movtrim(2)) = 0; 
    end
end

% Doing the actual dilation
masksdilated = imdilate(masksdilated > 0, strel('disk', p.npsize(2)));

% Pass
pass = ones(nmasks,1) > 0;

% Neuropil strel
npstrel = strel('disk', p.npsize(1));

% Loop through and fill in masks and neuropils
for i = 1 : nmasks
    % Get mask
    currentmask = masks == i;
        
    % Check edges
    if ~isempty(p.edges)
        currentmask(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2)) = 0;
    end
    
    % Check GRIN
    if p.GRIN
        currentmask = currentmask & p.grinface;
    end
    
    % Check nosoma
    if p.nosoma
        currentmask = currentmask & ~somaROIs;
    end
    
    % Dilate mask
    if p.binxy > 1
        currentmask = imresize(currentmask, p.binxy);
    end
    
    % Make sure that the mask dimension matches the original movie
    if trimmed
        if any(masktrim ~= 0)
            currentmask = currentmask(1:end-masktrim(1), 1:end-masktrim(2));
        end

        if movtrim(1) > 0
            currentmask(end+1 : end+movtrim(1), :) = 0; 
        end

        if movtrim(2) > 0
            currentmask(:, end+1 : end+movtrim(2)) = 0; 
        end
    end
    
    % Fill mask
    cellsort(i).mask = currentmask;
    
    % Fill neuropil (should be the right size)
    cellsort(i).neuropil = imdilate(currentmask, npstrel, 'same') & ~masksdilated;
    
    % Check edges
    if ~isempty(p.edges)
        cellsort(i).neuropil(1:p.edges(3), :) = 0;
        cellsort(i).neuropil(end-p.edges(4)+1:end, :) = 0;
        cellsort(i).neuropil(:, 1:p.edges(1)) = 0;
        cellsort(i).neuropil(:, end-p.edges(2)+1:end) = 0;
    end
    
    % Check GRIN
    if p.GRIN
        cellsort(i).neuropil = cellsort(i).neuropil & imresize(p.grinface, p.binxy);
    end
    
    % Check nosoma
    if p.nosoma
        cellsort(i).neuropil = cellsort(i).neuropil & imresize(~somaROIs, p.binxy);
    end
    
    % Properties
    properties = regionprops(currentmask, 'Area', 'Centroid', 'BoundingBox');
    if length(properties) == 1
        if properties.Area >= p.minarea
            cellsort(i).area = properties.Area;
            cellsort(i).centroid = round(properties.Centroid);
            cellsort(i).bbox = round(properties.BoundingBox);
        
        else
            pass(i) = false;
        end
    elseif p.allowsplit
        cellsort(i).area = [properties.Area];
        cellsort(i).centroid = round([properties.Centroid]);
        cellsort(i).bbox = round([properties.BoundingBox]);
        
        if isempty(cellsort(i).area) || sum(cellsort(i).area) < p.minarea
            pass(i) = false;
        end
    else
        pass(i) = false;
    end
end

% Remove cells that are in the wrong place
cellsort = cellsort(pass);

% Display
t = round(toc);
fprintf(' Done. Elapsed time = %i seconds.\n', t);

%% Pulling signals
% Signal
fprintf('Pulling raw signal...');
tic
cellsort = TiffPullSignalsCore(mov, cellsort, p.binxy);
t = round(toc);
fprintf(' Done. Elapsed time = %i seconds.\n', t);

% Remove those ROI with signal in neuropil that matches signal in ROI
fprintf('Processing correlated signal...');
tic
cellsort = sbxSignalsNeuropilCorrelation(cellsort);
t = round(toc);
fprintf(' Done. Elapsed time = %i seconds.\n', t);

% Get median-subtracted DFF Traces
if ~isempty(cellsort)
    fprintf('Calculating DFF...');
    tic
    cellsort = TiffSignalsWindowedDFF(cellsort, p.freq, p.dffwindow, p.percentile);
    t = round(toc);
    fprintf(' Done. Elapsed time = %i seconds.\n', t);   
end


%% Output
% Signal
sigstruct = struct('cellsort', cellsort, 'psig', p);

% Save file
if ~sigalt
    % Not alt signal
    save(sigpath, '-struct', 'sigstruct', '-v7.3');
else
    save(sigpath_alt, '-struct', 'sigstruct', '-v7.3');
end

% Output grinface
grinface = p.grinface;
end

