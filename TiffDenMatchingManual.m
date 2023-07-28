function TiffDenMatchingManual(mouse, date, runnum, optotune, varargin)
% TiffROIMatchingManual matches ROIs between experiments, using the first
% experiment as the main template. This function requires manual clicking.

%% Parse inputs
p = inputParser;

if nargin < 5
    varargin = {};
end

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'pmt', 1);
addOptional(p, 'sigsuffix', 'den');
addOptional(p, 'force', false);

% Bin
addOptional(p, 'binxy', 1); % Bin used for generating the to-seg image before

% Figure variables
addOptional(p, 'figpos', []);
addOptional(p, 'redalpha', 0.5);

% id field
addOptional(p, 'idfield', 'somaid');

% One to one pairing
addOptional(p, 'onetoone', false);

% String identifier (so other functions know which experiments are to be
% grouped together
addOptional(p, 'str', ''); % leave empty to copy str from main sig file
addOptional(p, 'justchangestr', false);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Figure locations
if isempty(p.figpos)
    figpos = [0 50 1300 850];
else
    figpos = p.figpos;
end

if ~ischar(date)
    date = num2str(date);
end

%% IO
% Signal file
if ~isempty(optotune) && ~isnan(optotune)
    sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
else
    sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
    [fp, ~, ~] = fileparts(sigpath);
    sigpath = fullfile(fp, sprintf('%s_%s_%03d.signals', mouse, date, runnum));
end
sigstruct = load(sigpath, '-mat', 'cellsort', 'fovstr');
ncells = length(sigstruct.cellsort);
nROIs_left = ncells;

% Make a file of all ROIs
ROIs = zeros(size(sigstruct.cellsort(1).mask));
for j = 1 : ncells
    ROIs(sigstruct.cellsort(j).mask) = j;
end
if p.binxy > 1
    ROIs = binxy(ROIs, p.binxy);
end
ROIs_left = ROIs;

% background file
[fp, ~, ~] = fileparts(sigpath);
if isempty(optotune) || isnan(optotune)
    meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runnum));
else
    meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runnum, optotune));
end
background = imread(meanpath);
    
% Den signal file
[fp_temp, fn_temp, ext_temp] = fileparts(sigpath);
fn_alt = sprintf('%s_%s%s', fn_temp, p.sigsuffix, ext_temp);
sigpath_den = fullfile(fp_temp, fn_alt);    
sigstruct_den = load(sigpath_den, '-mat', 'cellsort');
ncells_den = length(sigstruct_den.cellsort);
nROIs_left_den = ncells_den;

% No ROI
if ncells_den == 0
    warning('%s_%s_run%i cellsort is empty\n', mouse, date, runnum);
    return;
end

% Already done
if isfield(sigstruct_den(1).cellsort, p.idfield) && ~p.force
    redo = input('Matching already done. Redo (1 = yes, 0 = no): ');
    if redo ~= 1
        fprintf('%s_%s_run%i skipped\n', mouse, date, runnum);
        return;
    end
end

% Make a file of all ROIs den
ROIs_den = zeros(size(sigstruct_den.cellsort(1).mask));
for j = 1 : ncells_den
    ROIs_den(sigstruct_den.cellsort(j).mask) = j;
end
if p.binxy > 1
    ROIs_den = binxy(ROIs_den, p.binxy);
end
ROIs_left_den = ROIs_den;

% RGB den
RGB_den = repmat(mat2gray(background), [1 1 3]);
RGB_den(:,:,1) = (ROIs_den > 0) * p.redalpha;
RGB_den(:,:,3) = (ROIs > 0) * p.redalpha;

%% Intial plot
% ROI cell
matchingmat = -ones(ncells_den, 3);
matchingmat(:,1) = 1 : ncells_den;

% Figure den
hfig_den = figure('Position', figpos(1,:), 'name', sprintf('Run %i', runnum));
hrgb_den = imagesc(RGB_den);
set(hfig_den, 'MenuBar', 'none');
set(hfig_den, 'ToolBar', 'none');

% Current cell ID
cell_curr = 0;

% Getting user input
flagdone = false;

% Just changing strings
if p.justchangestr
    flagdone = true;
end

while ~flagdone
    % Advance cell count
    cell_curr = cell_curr + 1;
    
    % == Den ==
    figure(hfig_den);
    title(sprintf('%s: choose dendrite', date));
    
    % Figure out which ROI
    chosen_den = 0;
    while chosen_den == 0
        if nROIs_left_den > 1
            coord = round(ginput(1));
            if any(coord < 0) || coord(1)>hrgb_den.XData(2) || coord(2)>hrgb_den.YData(2)
                chosen_den = -1;
            else
                chosen_den = ROIs_left_den(coord(2), coord(1));
            end
        else
            chosen_den = unique(ROIs_left_den(ROIs_left_den > 0));
        end
    end
    matchingmat(cell_curr, 2) = chosen_den;
    hrgb_den.CData(:,:,1) = (ROIs_left_den == chosen_den) * p.redalpha;
    
    % == Main ==
    figure(hfig_den);
    title(sprintf('%s: choose soma', date));
    
    % Figure out which ROI
    chosen = 0;
    while chosen == 0
        coord = round(ginput(1));
        if any(coord < 0) || coord(1)>hrgb_den.XData(2) || coord(2)>hrgb_den.YData(2)
            chosen = -1;
        else
            chosen = ROIs_left(coord(2), coord(1));
        end
    end
    matchingmat(cell_curr, 3) = chosen;
    hrgb_den.CData(:,:,3) = (ROIs_left == chosen) * p.redalpha;
    
    % == Den ==
    % Remove ROI from pools
    if chosen_den > 0
        % Decraese ROI count by 1
        nROIs_left_den = nROIs_left_den - 1;

        % Remove ROI
        ROIs_left_den(ROIs_left_den == chosen_den) = 0;

        % Update RGB
        RGB_den(:,:,1) = (ROIs_left_den > 0) * p.redalpha;
    end

    % Update figure
%     hrgb_den.CData = RGB_den;
    
    % == Main ==
    % Remove ROI from pools
    if p.onetoone && chosen > 0
        % Decraese ROI count by 1
        nROIs_left = nROIs_left - 1;

        % Remove ROI
        ROIs_left(ROIs_left == chosen) = 0;

        % Update RGB
        RGB_den(:,:,1) = (ROIs_left > 0) * p.redalpha;
    end

    % Update figure
    hrgb_den.CData = RGB_den;

    % If no ROI left, it's done
    if nROIs_left == 0 || nROIs_left_den == 0
        flagdone = true;
    end
    
    % if both inputs are -1, it's done
    if chosen == -1 && chosen_den == -1
        flagdone = true;
    elseif chosen_den == -1
        % No dendrite selected but soma is selected, backtrack
        cell_curr = cell_curr - 1;
    end
    
    
end


%% Cleaning
% Matrix that denotes when ROIs are dropped
% Fix fitst column
if ~p.justchangestr
    matchingmat = matchingmat(1:cell_curr, :);
    matchingmat(:,1) = 1 : cell_curr;
end
close(hfig_den);

% Putting ids in dendrite struct
for i = 1 : ncells_den
    % Determine id
    switch p.idfield
        case 'somaid'
            id = matchingmat((matchingmat(:,2) == i), 3);
            if isempty(id)
                id = -1;
            end
    end
    
    % Put in id
    sigstruct_den.cellsort(i).(p.idfield) = id;
end

%% Saving
% Load up
loaded = load(sigpath_den, '-mat');
if ~p.justchangestr
    loaded.cellsort = sigstruct_den.cellsort;
    loaded.soma_ROIs = matchingmat;
end
if isempty(p.str)
    loaded.fovstr = sigstruct.fovstr;
else
    loaded.fovstr = p.str;
end
    
% Save
save(sigpath_den, '-struct', 'loaded', '-v7.3');

fprintf('%s_%s_run%i_OT%i_Soma matching done.\n', mouse, date, runnum, optotune);


end