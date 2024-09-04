function TiffROIMatchingManualXExpt(mousecell1, datecell1, runs1, optotunes1, ...
    mousecell2, datecell2, runs2, optotunes2, varargin)
% TiffROIMatchingManual matches ROIs between entirely different types of experiments, 
% using the first experiment as the main template. This function requires manual clicking.
% The experiments are separated into Group 1 and Group 2.
%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'useoptotune1', false); % Whether using optotune or not - Group 1
addOptional(p, 'useoptotune2', false); % Whether using optotune or not - Group 2
addOptional(p, 'pmt', 1);

% Bin
addOptional(p, 'binxy', 1); % Bin used for generating the to-seg image before

% Figure variables
addOptional(p, 'figpos', []);
addOptional(p, 'redalpha', 0.5);

% String identifier (so other functions know which experiments are to be
% grouped together
addOptional(p, 'str', '');

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs
% Number of sets
nset1 = length(mousecell1);
nset2 = length(mousecell2);
nset = nset1 + nset2;

% Figure locations
if isempty(p.figpos)
    if nset <= 6
        figpos = [100 560 560 420; 100 50 560 420; 700 560 560 420; 700 50 560 420;...
            1300 560 560 420; 1300 50 560 420];
    elseif nset <= 12
        figpos = [100 750 400 300; 100 400 400 300; 100 50 400 300; 550 750 400 300;...
            550 400 400 300; 550 50 400 300; 1000 750 400 300; 1000 400 400 300;...
            1000 50 400 300; 1450 750 400 300; 1450 400 400 300; 1450 50 400 300];
    end
else
    figpos = p.figpos;
end

% Turn everything into cell
if ~iscell(mousecell1)
    mousecell1 = {mousecell1};
end
if ~iscell(mousecell2)
    mousecell2 = {mousecell2};
end

if ~iscell(datecell1)
    datecell1 = {datecell1};
end
if ~iscell(datecell2)
    datecell2 = {datecell2};
end

if iscell(runs1)
    runs1 = cell2mat(runs1);
end
if iscell(runs2)
    runs2 = cell2mat(runs2);
end 


%% IO
% Initialize
ROI_struct = struct('mouse', [], 'date', [], 'run',[], 'optotune',[],...
    'cellsort', [], 'background', [], 'ROIs', [], 'ncells', []);
ROI_struct = repmat(ROI_struct, [nset, 1]);

% A cell to hold all the current ROIs
ROIs_left = cell(nset,1);

% A vector to hold the current number of ROIs
nROIs_left = zeros(nset,1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell = cell(nset,1);

hwait = waitbar(0, 'Preprocessing');
for i = 1 : nset
    % Update waitbar
    waitbar(i/nset, hwait, sprintf('Preprocessing %i/%i', i, nset));
    
    % Basic values
    if i <= nset1
        mouse = mousecell1{i};
        date = datecell1{i};
        runnum = runs1(i);
        optotune = optotunes1(i);

        % Load
        ROI_struct(i).mouse = mouse;
        ROI_struct(i).date = date;
        ROI_struct(i).run = runnum;
        ROI_struct(i).optotune = optotune;

        % Signal file
        if p.useoptotune1
            sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
        else
            sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
        end
        loaded = load(sigpath, '-mat', 'cellsort');
    else
        mouse = mousecell2{i-nset1};
        date = datecell2{i-nset1};
        runnum = runs2(i-nset1);
        optotune = optotunes2(i-nset1);

        % Load
        ROI_struct(i).mouse = mouse;
        ROI_struct(i).date = date;
        ROI_struct(i).run = runnum;
        ROI_struct(i).optotune = optotune;

        % Signal file
        if p.useoptotune2
            sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
        else
            sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
        end
        loaded = load(sigpath, '-mat', 'cellsort');
    end
    
    % cell sort
    ROI_struct(i).cellsort = loaded.cellsort;
    ncells = length(loaded.cellsort);
    ROI_struct(i).ncells = ncells;
    nROIs_left(i) = ncells;
    
    % Make a file of all ROIs (using xROI id)
    ROIs = zeros(size(loaded.cellsort(1).mask));
    for j = 1 : ncells
        ROIs(loaded.cellsort(j).mask) = loaded.cellsort(j).xrun_id;
    end
    if p.binxy > 1
        ROIs = binxy(ROIs, p.binxy);
    end
    ROI_struct(i).ROIs = ROIs;
    ROIs_left{i} = ROIs;
    
    % background file
    [fp, ~, ~] = fileparts(sigpath);
    if i <= nset1
        if ~p.useoptotune1
            meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runnum));
        else
            meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runnum, optotune));
        end
    else
        if ~p.useoptotune2
            meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runnum));
        else
            meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runnum, optotune));
        end
    end
    background = imread(meanpath);
    ROI_struct(i).background = background;
    
    % RGB
    RGB = repmat(mat2gray(background), [1 1 3]);
    RGB(:,:,1) = (ROIs > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell{i} = RGB;
end
close(hwait)

%% Intial plot
% Cells to get all the handles
hfigs = cell(nset,1);
hrgb = cell(nset,1);

% ROI cell
matchingmat = ones(ROI_struct(1).ncells, nset+1);
matchingmat(:,1) = 1 : ROI_struct(1).ncells;

for i = 1 : nset
    % Figure
    if i <= nset1
        hfigs{i} = figure('Position', figpos(i,:), 'name', 'Expt 1');
    else
        hfigs{i} = figure('Position', figpos(i,:), 'name', 'Expt 2');
    end
    hrgb{i} = imagesc(RGB_cell{i});
    if i <= nset1
        title(sprintf('Date %s', datecell1{i}));
    else
        title(sprintf('Date %s', datecell2{i-nset1}));
    end
    
    set(hfigs{i}, 'MenuBar', 'none');
    set(hfigs{i}, 'ToolBar', 'none');
end

% Current cell ID
cell_curr = 0;

% Getting user input
flagdone = false;
while ~flagdone
    % Initialize chosen vector
    chosen = zeros(nset,1);
    
    % Advance cell count
    cell_curr = cell_curr + 1;
    
    % Pick 1 time per experiment
    Expt1_picked = false;
    Expt2_picked = false;
    
    % Get inputs
    for i = 1 : nset
        % Get input
        figure(hfigs{i});
        
        % Skip
        skip = false;
        
        % Auto skip if none left
        if nROIs_left(i) == 0
            skip = true;
            chosen(i) = -1;
        end
        
        % Auto skip if this experiment is already done
        if i <= nset1 && Expt1_picked
            skip = true;
            chosen(i) = chosen(i-1);
        end
        if i > nset1 && Expt2_picked
            skip = true;
            chosen(i) = chosen(i-1);
        end
        
        % Figure which one was chosen
        while chosen(i) == 0
            % Figure out which ROI
            coord = round(ginput(1));
            
            if any(coord < 0) || coord(1)>hrgb{i}.XData(2) || coord(2)>hrgb{i}.YData(2)
                skip = true;
                chosen(i) = -1;
            else
                chosen(i) = ROIs_left{i}(coord(2), coord(1));
                
                % Tell subsequent runs that this experiment is already
                % picked
                if i <= nset1
                    Expt1_picked = true;
                else
                    Expt2_picked = true;
                end
            end
        end

        % Update cdata
        if skip && nROIs_left(i) == 0
             hrgb{i}.CData(:,:,1) = 0;
        else
             hrgb{i}.CData(:,:,1) = (ROIs_left{i} == chosen(i)) * p.redalpha;
        end
       
    end

    % Apply inputs
    for i = 1 : nset
        % Apply cellsort
        if chosen(i) > 0
            % Find what is the current id (convert from xrun id)
            k = [ROI_struct(i).cellsort(:).xrun_id] == chosen(i);
            if any(k)
                % This xrun id exists for this run
                ROI_struct(i).cellsort(k).xexpt_id = cell_curr;
            else
                % This xrun id does not exist for this run
                chosen(i) = -1;
            end
        end
        
        % Update cell
        matchingmat(cell_curr, i+1)= chosen(i);
    end

    % Remove ROI from pools
    for i = 1 : nset
        if chosen(i) > 0
            % Decraese ROI count by 1
            nROIs_left(i) = nROIs_left(i) - 1;

            % Remove ROI
            ROIs_left{i}(ROIs_left{i} == chosen(i)) = 0;

            % Update RGB
            RGB_cell{i}(:,:,1) = (ROIs_left{i} > 0) * p.redalpha;
        end
        
        % Update figure
        hrgb{i}.CData = RGB_cell{i};
    end

    % If no ROI left, it's done
    if sum(nROIs_left(1:nset1)) == 0
        flagdone = true;
    end
    if sum(nROIs_left(nset1+1 : end)) == 0
        flagdone = true;
    end
end


%% Cleaning
% Matrix that denotes when ROIs are dropped
% Fix fitst column
matchingmat(:,1) = 1 : cell_curr;

for i = 1 : nset
    close(hfigs{i});
end

%% Saving
hwait = waitbar(0, 'Saving');
for i = 1 : nset
    % Update waitbar
    waitbar(i/nset, hwait, sprintf('Saving %i/%i', i, nset));
    
    if i <= nset1
        % Basic values
        mouse = mousecell1{i};
        date = datecell1{i};
        runnum = runs1(i);
        optotune = optotunes1(i);

        % Signal file
        if p.useoptotune1
            sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
        else
            sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
        end
    else
        % Basic values
        mouse = mousecell2{i-nset1};
        date = datecell2{i-nset1};
        runnum = runs2(i-nset1);
        optotune = optotunes2(i-nset1);

        % Signal file
        if p.useoptotune2
            sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
        else
            sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
        end
    end
    
    % Update structures
    loaded = load(sigpath, '-mat');
    loaded.cellsort = ROI_struct(i).cellsort;
    loaded.xrun_ROIs = matchingmat;
    loaded.fovstr = p.str;
    
    % Save
    save(sigpath, '-struct', 'loaded', '-v7.3');
end
close(hwait)
disp('xexpt matching done.')

end