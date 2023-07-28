function TiffROIMatchingManual(mousecell, datecell, runs, optotunes, varargin)
% TiffROIMatchingManual matches ROIs between experiments, using the first
% experiment as the main template. This function requires manual clicking.

%% Parse inputs
p = inputParser;

% Path variables
addOptional(p, 'server', 'nasquatch');
addOptional(p, 'useoptotune', false); % Whether using optotune or not
addOptional(p, 'pmt', 1);

% Bin
addOptional(p, 'binxy', 1); % Bin used for generating the to-seg image before

% Figure variables
addOptional(p, 'figpos', []);
addOptional(p, 'redalpha', 0.5);

% String identifier (so other functions know which experiments are to be
% grouped together
addOptional(p, 'str', '');
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
    if length(mousecell) <= 6
        figpos = [100 560 560 420; 100 50 560 420; 700 560 560 420; 700 50 560 420;...
            1300 560 560 420; 1300 50 560 420];
    elseif length(mousecell) <= 12
        figpos = [100 750 400 300; 100 400 400 300; 100 50 400 300; 550 750 400 300;...
            550 400 400 300; 550 50 400 300; 1000 750 400 300; 1000 400 400 300;...
            1000 50 400 300; 1450 750 400 300; 1450 400 400 300; 1450 50 400 300];
    elseif length(mousecell) <= 63
        % Auto determination
        screenheight = 1080;
        screenwidth = 1920;
        
        if length(mousecell) <= 20
            width = 356; %356 (20) | 285 (30) | 237 (48)
            height = 267; %267 (20) | 214 (30) | 178 (48)
        elseif length(mousecell) <= 30 % 5 * 6
            width = 285; %356 (20) | 285 (30) | 237 (48)
            height = 214; %267 (20) | 214 (30) | 178 (48)
        elseif length(mousecell) <= 48 % 6 * 8
            width = 285; %356 (20) | 285 (30) | 237 (48)
            height = 214; %267 (20) | 214 (30) | 178 (48)
        elseif length(mousecell) <= 63 % 7 x 9
            width = 200; %356 (20) | 285 (30) | 237 (48)
            height = 150; %267 (20) | 214 (30) | 178 (48)
        end
        
        % Panels to fit
        npanelsy = floor((screenheight-10) / height);
        npanelsx = floor((screenwidth-10) / width);
        
        % Give a bit more gap if possible
        heightgap = (screenheight-10) / npanelsy;
        widthgap = (screenwidth-10) / npanelsx;
        
        % Get the positions
        x = 10 : widthgap : screenwidth-width;
        y = screenheight-height : -heightgap : 10;
        xvec = ones(length(y),1) * x;
        yvec = y' * ones(1, length(x));
        
        figpos = zeros(length(x)*length(y), 4);
        figpos(:,1) = xvec(:);
        figpos(:,2) = yvec(:);
        figpos(:,3) = width;
        figpos(:,4) = height;
    end
else
    figpos = p.figpos;
end

% Optotune number to vector
if length(optotunes) == 1
    optotunes = ones(length(mousecell),1) * optotunes;
end

% Turn everything into cell
if ~iscell(mousecell)
    mousecell = {mousecell};
end

if ~iscell(datecell)
    datecell = {datecell};
end

if iscell(runs)
    runs = cell2mat(runs);
end
    
% Number of sets
nset = length(mousecell);

%% IO
% Initialize
ROI_struct = struct('mouse', [], 'date', [], 'run',[], 'optotune',[],...
    'cellsort', [], 'background', [], 'ROIs', [], 'ncells', []);
ROI_struct = repmat(ROI_struct, [nset, 1]);

% A cell to hold all the current ROIs
ROIs_left = cell(length(runs),1);

% A vector to hold the current number of ROIs
nROIs_left = zeros(length(runs),1);

% A cell to hold all the current RGBs (R = ROI, G = photon image)
RGB_cell = cell(length(runs),1);

hwait = waitbar(0, 'Loading');
for i = 1 : nset
    % waitbar
    waitbar(i/nset, hwait, sprintf('Loading %i/%i.', i, nset));
    
    % Basic values
    mouse = mousecell{i};
    date = datecell{i};
    if ~ischar(date)
        date = num2str(date);
    end
    runnum = runs(i);
    optotune = optotunes(i);
    
    % Load
    ROI_struct(i).mouse = mouse;
    ROI_struct(i).date = date;
    ROI_struct(i).run = runnum;
    ROI_struct(i).optotune = optotune;
    
    % Signal file
    if p.useoptotune && ~isempty(optotune) && ~isnan(optotune)
        sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
    else
        sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
    end
    loaded = load(sigpath, '-mat', 'cellsort');
    
    % cell sort
    ROI_struct(i).cellsort = loaded.cellsort;
    ncells = length(loaded.cellsort);
    ROI_struct(i).ncells = ncells;
    nROIs_left(i) = ncells;
    
    % Make a file of all ROIs
    ROIs = zeros(size(loaded.cellsort(1).mask));
    for j = 1 : ncells
        ROIs(loaded.cellsort(j).mask) = j;
    end
    if p.binxy > 1
        ROIs = binxy(ROIs, p.binxy);
    end
    ROI_struct(i).ROIs = ROIs;
    ROIs_left{i} = ROIs;
    
    % background file
    [fp, ~, ~] = fileparts(sigpath);
    if ~p.useoptotune || isempty(optotune) || isnan(optotune)
        meanpath = fullfile(fp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, runnum));
    else
        meanpath = fullfile(fp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, runnum, optotune));
    end
    background = imread(meanpath);
    ROI_struct(i).background = background;
    
    % RGB
    RGB = repmat(mat2gray(background), [1 1 3]);
    RGB(:,:,1) = (ROIs > 0) * p.redalpha;
    RGB(:,:,3) = 0;
    RGB_cell{i} = RGB;
end
close(hwait);

%% Intial plot
% Cells to get all the handles
hfigs = cell(nset,1);
hrgb = cell(nset,1);

% ROI cell
matchingmat = ones(ROI_struct(1).ncells, nset+1);
matchingmat(:,1) = 1 : ROI_struct(1).ncells;

for i = 1 : length(runs)
    % Figure
    hfigs{i} = figure('Position', figpos(i,:), 'name', sprintf('Run %i', runs(i)));
    hrgb{i} = imagesc(RGB_cell{i});
    date = datecell{i};
    if ~ischar(date)
        date = num2str(date);
    end
    title(sprintf('Date %s', date));
    
    set(hfigs{i}, 'MenuBar', 'none');
    set(hfigs{i}, 'ToolBar', 'none');
end

% Current cell ID
cell_curr = 0;

% Getting user input
flagdone = false;

% Just changing strings
if p.justchangestr
    flagdone = true;
end

while ~flagdone
    % Initialize chosen vector
    chosen = zeros(nset,1);
    
    % Advance cell count
    cell_curr = cell_curr + 1;
    
    % If only 1 panel left, autochoose
    if sum(nROIs_left > 0) == 1
        panelleft = find(nROIs_left > 0);
        chosen(panelleft) = max(max(ROIs_left{panelleft}));
        autochoose = true;
    else
        autochoose = false;
    end
    
    % Get inputs
    for i = 1 : length(runs)
        % Get input
        figure(hfigs{i});
        
        % Skip
        skip = false;
        
        % Auto skip if none left
        if nROIs_left(i) == 0
            skip = true;
            chosen(i) = -1;
        end
        
        % Autochoose because only 1 panel left
        if autochoose
            skip = true;
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
            end
        end

        % Update cdata
        if skip
             hrgb{i}.CData(:,:,1) = 0;
        else
             hrgb{i}.CData(:,:,1) = (ROIs_left{i} == chosen(i)) * p.redalpha;
        end
       
    end

    % Apply inputs
    for i = 1 : length(runs)
        % Apply cellsort
        if chosen(i) > 0
            ROI_struct(i).cellsort(chosen(i)).xrun_id = cell_curr;
        end
        
        % Update cell
        matchingmat(cell_curr, i+1)= chosen(i);
    end

    % Remove ROI from pools
    for i = 1 : length(runs)
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
    if sum(nROIs_left) == 0
        flagdone = true;
    end
end


%% Cleaning
% Matrix that denotes when ROIs are dropped
% Fix fitst column
if ~p.justchangestr
    matchingmat(:,1) = 1 : cell_curr;
end
for i = 1 : nset
    close(hfigs{i});
end

%% Saving
hwait = waitbar(0, 'Saving');
for i = 1 : nset
    % Save
    waitbar(i/nset, hwait, sprintf('Saving %i/%i.', i, nset));
    
    % Basic values
    mouse = mousecell{i};
    date = datecell{i};
    runnum = runs(i);
    optotune = optotunes(i);
    
    % Signal file
    if p.useoptotune && ~isempty(optotune) && ~isnan(optotune)
        sigpath = sbxPath(mouse, date, runnum, 'OTsig', 'server', p.server, 'optotune', optotune);
    else
        sigpath = sbxPath(mouse, date, runnum, 'signals', 'server', p.server);
    end
    
    % Update structures
    loaded = load(sigpath, '-mat');
    if ~p.justchangestr
        loaded.cellsort = ROI_struct(i).cellsort;
        loaded.xrun_ROIs = matchingmat;
    end
    loaded.fovstr = p.str;
    
    % Save
    save(sigpath, '-struct', 'loaded', '-v7.3');
end
disp('xrun matching done.')
close(hwait)

end