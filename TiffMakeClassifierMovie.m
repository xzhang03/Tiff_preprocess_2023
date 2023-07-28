function TiffMakeClassifierMovie(mousecell, datecell, runvec, optotunecell, varargin)
%TiffMakeClassifierMovie makes movie based on classifier result
% TiffMakeClassifierMovie(mousecell, datecell, runvec, optotunecell, varargin)

%% Parse inputs
p = inputParser;
    
% Lab variables variables
addOptional(p, 'server', []);  % Add in the server name as a string
addOptional(p, 'force', false); % Force reclassify or not
addOptional(p, 'fps', []); % Fps - must supply
addOptional(p, 'pmt', 1);
addOptional(p, 'useoptotune', false);
addOptional(p, 'movtype', 'tiff_demonsreg');

% Which value in pass to use
addOptional(p, 'passval', 1);
addOptional(p, 'failval', []); % Leave empty to make it anything else that doesn't pass

% Motion correlation
addOptional(p, 'usemocorr', true);
addOptional(p, 'mocorrthresh', 0.5); % In abs value

% Neuropil correlation
addOptional(p, 'usenpcorr', false);
addOptional(p, 'npcorrthresh', 0.7); % In signed value

% Preprocessing
addOptional(p, 'smooth', false); % Smoothing
addOptional(p, 'smoothwin', 3);
addOptional(p, 'movmed', false); % Moving median
addOptional(p, 'movmedwin', 3);

% Movie parameters
addOptional(p, 'sizevec', [20 20]); % rows and columns from centroid. E.g., 10 means -10 to +10 pixels.
addOptional(p, 'prestim', 20); % Pre-stim window
addOptional(p, 'poststim', 50); % Post-stim dinwo
addOptional(p, 'usesigonsets', false); % Use signal onsets as opposed to event onsets

% dff
addOptional(p, 'dffmovie', true); % dff movie
addOptional(p, 'dfftime', []); % In seconds, leave empty to make it the same as the prestim window
addOptional(p, 'dffbeforeaverage', false); % Dff before averaging
addOptional(p, 'znorm', false); % z normalize

% xrun parameters
addOptional(p, 'usexrun', true); % Use xrun ROI ids
addOptional(p, 'fovs', 1); % Number of fovs, used for initialization only. Doesn't have to be accurate

% Binning parameters
addOptional(p, 'binxy', 1); 
addOptional(p, 'bint', 1);

% Movie output
addOptional(p, 'nfperrow', []);
addOptional(p, 'concatenatedir', 2); % 2 = horizontally, 1 = vertically
addOptional(p, 'tiling', 'dense'); % can be dense or matching

% Markers
addOptional(p, 'usewhitebox', true);
addOptional(p, 'whiteboxsize', 10);
addOptional(p, 'whiteboxvalue', [60000 3]); % First value is absolute in non-dff mode. Second value in diff mode.
addOptional(p, 'separationbar', true);
addOptional(p, 'separationbarsize', 10);

% Mask
addOptional(p, 'showmask', true); % A contour of mask
addOptional(p, 'maskvalue', [500, 3]); % First value is absolute in non-dff mode. Second value in diff mode.

% Output path
addOptional(p, 'outputpath', '');
addOptional(p, 'outputfn', '');

% Save parameters
addOptional(p, 'savepars', true);

% Save materials
addOptional(p, 'savematerials', false);

% Debug (only load 100 frames);
addOptional(p, 'debug', false);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

%% Clean up inputs and initializing
% Turn everything into cell
if ~iscell(mousecell)
    mousecell = {mousecell};
end

if ~iscell(datecell)
    datecell = {datecell};
end

if iscell(runvec)
    runvec = cell2mat(runvec);
end
    
if ~iscell(optotunecell)
    optotunecell = {optotunecell};
end

% Numbet of datasets
if ~p.useoptotune
    % No optotune situation
    nset = length(mousecell);
    
    % Optotune cell
    optotunecell = cell(nset,1);
else
    % With optotune, remake the cells
    nset = 0;
    
    % Remake (initialize with 100)
    mousecell2 = cell(100,1);
    datecell2 = cell(100,1);
    runvec2 = zeros(100,1);
    optotunecell2 = cell(100,1);
    
    for i = 1 : length(mousecell)
        % N
        nOT = length(optotunecell{i});
        
        % Populate
        mousecell2(nset+1 : nset+nOT) = mousecell(i);
        datecell2(nset+1 : nset+nOT) = datecell(i);
        runvec2(nset+1 : nset+nOT) = runvec(i);
        optotunecell2(nset+1 : nset+nOT) = num2cell(optotunecell{i});
        
        nset = nset + nOT;
    end
    
    % Return
    mousecell = mousecell2(1:nset);
    datecell = datecell2(1:nset);
    runvec = runvec2(1:nset);
    optotunecell = optotunecell2(1:nset);    
end

if isempty(p.outputpath) || isempty(p.outputfn)
    disp('Not enough output path information is specified')
    return
end

if exist(fullfile(p.outputpath, p.outputfn), 'file') && ~p.force
    disp('Output movie already exists. Exiting.')
    return
end

% Vector of cell counts
ncellvec = zeros(nset, 1);

% Vector of all trial counts
ntrialvec = zeros(nset, 1);

% Cell of passfails
passcell = cell(nset, 1);

% Cell of ons and offs
onoffcell = cell(nset,2);

% Cell of motion corrs
mocorrcell = cell(nset, 1);

% Cell of np corrs
npcorrcell = cell(nset, 1);

% Centroid cenn
cencell = cell(nset,1);

% Delays (using signal triggers)
delaycell = cell(nset, 1);

% Z
zcell = cell(nset, 1);

%% Process the signal files
% If cross run store strings to fovs
if p.usexrun
    % Strings
    xrunstrs = cell(p.fovs, 1);
else
    xrunstrs = 'Noxrun';
end

% Each run is considered a new fov in the non-xrun mode
% Counts
ncell_xrun = zeros(p.fovs, 1);
ifov = 0;

% fov vector
fovvec = zeros(nset,1);

% xrun cell id (initialize with 100)
cellid_xrun = zeros(100,nset);

% masks (initialize with 100)
if p.showmask
    masks_cell = cell(100, nset);
end

for i = 1 : nset
    % Display
    fprintf('Loading data set %i/%i.\n', i, nset);

    % Data path
    if p.useoptotune
        sigpath = sbxPath(mousecell{i}, datecell{i}, runvec(i), 'OTsig', 'server', p.server, 'optotune', optotunecell{i});
    else
        sigpath = sbxPath(mousecell{i}, datecell{i}, runvec(i), 'signals', 'server', p.server);
    end
    
    sigstruct = load(sigpath, '-mat');
    
    % Basic paramteres
    ncellvec(i) = size(sigstruct.trigtensor, 3);
    ntrialvec(i) = size(sigstruct.trigtensor, 2);
    onoffcell{i,1} = sigstruct.nons;
    onoffcell{i,2} = sigstruct.noffs;
    passcell{i} = reshape(sigstruct.Group, ntrialvec(i), ncellvec(i));
    mocorrcell{i} = sigstruct.mocorr;
    npcorrcell{i} = sigstruct.npcorr;
    cencell{i} = reshape([sigstruct.cellsort(:).centroid],2,[])';
    
    % Delay onsets (use signal triggers as opposed to event triggers)
    if p.usesigonsets && isfield(sigstruct, 'onsets')
        delaycell{i} = reshape(sigstruct.onsets, ntrialvec(i), ncellvec(i));
    end
    
    % Z normalize
    if p.znorm
        zcell{i} = squeeze(std(sigstruct.trigtensor,1));
    end
        
    % xruncell id
    if p.usexrun
        cellid_xrun(1:ncellvec(i),i) = [sigstruct.cellsort(:).xrun_id];
    else
        cellid_xrun(1:ncellvec(i),i) = 1 : ncellvec(i);
    end
    
    % Use the first run to initialize xrun
    if p.usexrun
        % Comparing fov strings
        cp = strcmp(sigstruct.fovstr, xrunstrs);
        
        if ~any(cp)
            % New fov
            ifov = ifov + 1;
            
            % Store cell counts and strings
            ncell_xrun(ifov) = size(sigstruct.xrun_ROIs,1);
            xrunstrs{ifov} = sigstruct.fovstr;
            
            % Vector
            fovvec(i) = ifov;
        else
            % Existing fov
            fovvec(i) = find(cp, 1);
        end
    else
        % If not xrun treat each run as a new fov
        ncell_xrun(i) = ncellvec(i);
        ifov = ifov + 1;
        fovvec(i) = ifov;
    end
    
    % masks
    if p.showmask
        movdim = size(sigstruct.cellsort(1).mask);
        
        for j = 1 : ncellvec(i)
            % Find where the cellis
            cen = cencell{i}(j,:);
            cenvec = [cen(2)-p.sizevec(1), cen(2)+p.sizevec(1), cen(1)-p.sizevec(2), cen(1)+p.sizevec(2)];

            % See if the box is out of range
            if cenvec(1) <= 0
                ys = 1 - cenvec(1);
                cenvec = cenvec + [ys ys 0 0];
            end
            if cenvec(2) > movdim(1)
                ys = movdim(1) - cenvec(2);
                cenvec = cenvec + [ys ys 0 0];
            end
            if cenvec(3) <= 0
                xs = 1 - cenvec(3);
                cenvec = cenvec + [0 0 xs xs];
            end
            if cenvec(4) > movdim(2)
                xs = movdim(2) - cenvec(4);
                cenvec = cenvec + [0 0 xs xs];
            end
            
            % Get mask
            mask = sigstruct.cellsort(j).mask(cenvec(1) : cenvec(2), cenvec(3) : cenvec(4));
            
            % Binxy if needed
            if p.binxy > 1
                mask = binxy(mask, p.binxy);
            end
            
            % Save
            masks_cell{j,i} = edge(mask);
        end
    end
end

% Total number of cells
% Remove excess
nfov = ifov;
ncell_xrun = ncell_xrun(1:nfov);
ncell = sum(ncell_xrun);
cellid_xrun = cellid_xrun(1:max(ncellvec),:);

if p.showmask
    masks_cell = masks_cell(1:max(ncellvec),:);
end

if isempty(p.nfperrow)
    p.nfperrow = ceil(sqrt(ncell));
end

% Display
fprintf('Found %i FOVs with a total of %i cells.\n', nfov, ncell);

%% Pass and fail
% Initialize
% Max ncells in a run
mcells = max(ncell_xrun);

% pass and fail cells
passindcell = cell(mcells, nset);
failindcell = cell(mcells, nset);

% number vectors
np_mat = zeros(mcells, nfov);
nf_mat = zeros(mcells, nfov);

% xrun masks
if p.showmask
    masks_xrun = cell(mcells, nfov);
end

% Summing
np_total = 0;
nf_total = 0;
throw_total = 0;

for i = 1 : nset
    % Which fov
    ifov = fovvec(i);
    
    % Get the right matrices
    pass = passcell{i};
    mocorr = mocorrcell{i};
    npcorr = npcorrcell{i};
    ntrials = ntrialvec(i);
    
    for j = 1 : ncellvec(i)
        % xrun cellid
        cellid_curr = cellid_xrun(j,i);
        
        % Cells to consider (passed mocorr test)
        if p.usemocorr
            c = abs(mocorr(:,j)) <= p.mocorrthresh;
        else
            c = ones(ntrials,1) > 0;
        end
        
        % Cells to consider based on neuropil signal correlation
        if p.usenpcorr
            c = (npcorr(:,j) <= p.npcorrthresh) & c;
        end

        % Cells in the categories
        passc = (pass(:,j) == p.passval) & c;
        if ~isempty(p.failval)
            if length(p.failval) == 2
                % Require for 2 consequtive trials
                failc = (pass(:,j) == p.failval(1)) & c;
                failc2 = (pass(2:end,j) == p.failval(2)) & c(2:end);
                failc = failc & [failc2; false];
            elseif length(p.failval) == 1
                failc = (pass(:,j) == p.failval) & c;
            end
        else
            failc = (pass(:,j) ~= p.passval) & c;
        end
        
        % Save indices
        passindcell{cellid_curr, i} = passc;
        failindcell{cellid_curr, i} = failc;
        
        % numbers
        np = sum(passc);
        nf = sum(failc);
        
        % Add to sum
        np_total = np_total + np;
        nf_total = nf_total + nf;
        throw_total = throw_total + sum(~c);
        
        % Save numbers
        np_mat(cellid_curr, ifov) = np_mat(cellid_curr, ifov) + np;
        nf_mat(cellid_curr, ifov) = nf_mat(cellid_curr, ifov) + nf;
        
        % Get mask
        if p.showmask
            % Get mask
            mask = masks_cell{j, i};
            
            % Put it in the right xrun place
            if isempty(masks_xrun{j, i})
                masks_xrun{j, i} = mask;
            else
                masks_xrun{j, i} = masks_xrun{j, i} & mask;
            end
        end
    end
end

% Display
fprintf('Found %i passed trials and %i fails. Pass rate = %0.2f.\n', np_total, nf_total, np_total/(np_total+nf_total));
fprintf('Found %i trials that did not pass motion check. Pass rate = %0.2f.\n', throw_total, 1 - throw_total/sum(ncellvec.*ntrialvec));

%% Initialize movie
% Grid (two grids: responsive - unreponsive)
grid = [ceil(ncell/p.nfperrow), p.nfperrow];

% Frames before stim
preframes = ceil(p.prestim * p.fps);

% Bin
preframes2 = floor(preframes / p.bint);

% dff frames
if isempty(p.dfftime) && p.dffmovie
    dffframes = preframes;
elseif p.dffmovie
    dffframes = ceil(p.dfftime * p.fps);
end

if p. dffmovie
    % Bin
    dffframes2 = floor(dffframes / p.bint);
end

% Frames after stim
postframes = ceil(p.poststim * p.fps);

% Length
l = preframes + postframes + 1;

% Bin
l2 = floor(l / p.bint);

% Size
panelsize = p.sizevec * 2 + 1;

% Bin
panelsize2 = floor(panelsize / p.binxy);

%% Initialize cells to contain all the triggered movie fragments across runs
% Cell
mov4d_p_cell = cell(size(np_mat));
mov4d_f_cell = cell(size(nf_mat));

% Initialize 4D matrices to save triggered movies
for i = 1 : nfov
    for j = 1 : mcells
        if np_mat(j, i) > 0
            mov4d_p_cell{j, i} = zeros(panelsize2(1), panelsize2(2), l2, np_mat(j, i));
        end
        if nf_mat(j, i) > 0
            mov4d_f_cell{j, i} = zeros(panelsize2(1), panelsize2(2), l2, nf_mat(j, i));
        end
    end
end

% Initialize a matrix to keep tracl of frag inds
ip_mat = zeros(size(np_mat));
if_mat = zeros(size(nf_mat));

%% Load movie
for i = 1 : nset
    % fov
    ifov = fovvec(i);
    
    % Values
    mouse = mousecell{i};
    date = datecell{i};
    runnum = runvec(i);
    optotune = optotunecell{i};
    ncells = ncellvec(i);
%     ntrials = ntrialvec(i);
    nons = onoffcell{i,1};
    noffs = onoffcell{i,2};
    boxlength = round(mean(noffs - nons) / p.bint);
    cens = cencell{i};
    
    % Delay onset
    if p.usesigonsets
        delays = delaycell{i};
    end
    
    % Znorm
    if p.znorm
        zmat = zcell{i};
    end
    
    % Display
    fprintf('Processing Movie %i/%i, FOV %i, %i cells.\n', i, nset, ifov, ncells)
    
    % Debug mode
    if p.debug
        boxlength = 4;
        postframes = 10;
        preframes = 10;
        preframes2 = floor(preframes/p.bint);
        l = preframes + postframes + 1;
        l2 = floor(l / p.bint);
        p.smoothwin = 1;
        p.movmedwin = 1;
    end
    
    % Path
    if p.useoptotune
        tiffpath = sbxPath(mouse, date, runnum, p.movtype, 'server', p.server, 'pmt', p.pmt,...
            'optotune', optotune);
    else
        tiffpath = sbxPath(mouse, date, runnum, p.movtype, 'server', p.server, 'pmt', p.pmt,...
            'optotune', []);
    end

    % Load movie
    mov = readtiff(tiffpath);
    movdim = size(mov);
%     mov = mean(mov,3);
    
    % Loop through cells
    for j = 1 : ncells
        % xrun cellid
        cellid_curr = cellid_xrun(j,i);
        
        % nons for this cell
        passc = passindcell{cellid_curr, i};
        failc = failindcell{cellid_curr, i};
        if p.usesigonsets
            nons2 = nons + delays(:,j);
            nonsp = nons2(passc);
            nonsf = nons(failc);
        else
            nonsp = nons(passc);
            nonsf = nons(failc);
        end
        
        % znorm
        if p.znorm
            zmp = zmat(passc,j);
            zmp_m = median(zmp);
            zmf = zmat(failc,j);
            zmf_m = median(zmf);
        end
        
        % Pass and fail trials of the current movie
        np = sum(passc);
        nf = sum(failc);
        
        % indmat for this cell
        if p.smooth || p.movmed
            % Trim
            trim = max(p.smoothwin * p.smooth, p.movmedwin * p.movmed);
            trim = ceil(trim * p.fps);
            
            % Indmats
            % Ind mat for pass
            lplus = l + trim * 2;
            indmatp1 = nonsp * ones(1,lplus);
            indmatp2 = ones(np,1) * (-preframes-trim : postframes+trim);
            indmatp = indmatp1 + indmatp2;
            
            % Ind mat for fail
            indmatf1 = nonsf * ones(1,l + trim * 2);
            indmatf2 = ones(nf,1) * (-preframes-trim : postframes+trim);
            indmatf = indmatf1 + indmatf2;
        else
            % Ind mat for pass
            indmatp1 = nonsp * ones(1,l);
            indmatp2 = ones(np,1) * (-preframes : postframes);
            indmatp = indmatp1 + indmatp2;
            
            % Ind mat for fail
            indmatf1 = nonsf * ones(1,l);
            indmatf2 = ones(nf,1) * (-preframes : postframes);
            indmatf = indmatf1 + indmatf2;
        end
        
        % Find where the cellis
        cen = cens(j,:);
        cenvec = [cen(2)-p.sizevec(1), cen(2)+p.sizevec(1), cen(1)-p.sizevec(2), cen(1)+p.sizevec(2)];
        
        % See if the box is out of range
        if cenvec(1) <= 0
            ys = 1 - cenvec(1);
            cenvec = cenvec + [ys ys 0 0];
        end
        if cenvec(2) > movdim(1)
            ys = movdim(1) - cenvec(2);
            cenvec = cenvec + [ys ys 0 0];
        end
        if cenvec(3) <= 0
            xs = 1 - cenvec(3);
            cenvec = cenvec + [0 0 xs xs];
        end
        if cenvec(4) > movdim(2)
            xs = movdim(2) - cenvec(4);
            cenvec = cenvec + [0 0 xs xs];
        end
        
        % Loop through pass trials
        if np > 0
            % Get initial frag ind
            fi = ip_mat(cellid_curr, ifov);
            
            for k = 1 : np
                % Advance frag id
                fi = fi + 1;
                
                % Get a slice
                movslice = mov(cenvec(1) : cenvec(2), cenvec(3) : cenvec(4), indmatp(k,:));

                % Binxy if needed
                if p.binxy > 1
                    movslice = binxy(movslice, p.binxy);
                end

                % Smooth if needed
                if p.smooth
                    movslice = movmean(movslice, p.smoothwin * p.fps, 3);
                end

                % Movmed if needed
                if p.movmed
                    movslice = movmedian(movslice, p.movmedwin * p.fps, 3);
                end

                % Remove trim
                if p.smooth || p.movmed
                    movslice = movslice(:,:,trim+1 : end-trim);
                end

                % Bint if needed
                if p.bint > 1
                    movslice = bint(movslice, p.bint);
                end
                
                if p.dffbeforeaverage && p.dffmovie
                    % Average frame
                    movslice_base = mean(movslice(:,:,1:dffframes2), 3);
                    movslice = movslice ./ repmat(movslice_base, [1 1 l2]) - 1;
                end
                
                if p.znorm
                    movslice = movslice / zmp(k) * zmp_m;
                end
                
                % Save
                mov4d_p_cell{cellid_curr, ifov}(:,:,:,fi) = movslice;
            end
            
            % Save frag ind
            ip_mat(cellid_curr, ifov) = fi;
        end
        
        % Loop through fail trials
        if nf > 0
            % Get initial frag ind
            fi = if_mat(cellid_curr, ifov);
            
            for k = 1 : nf
                % Advance frag id
                fi = fi + 1;
                
                % Get a slice
                movslice = mov(cenvec(1) : cenvec(2), cenvec(3) : cenvec(4), indmatf(k,:));

                % Binxy if needed
                if p.binxy > 1
                    movslice = binxy(movslice, p.binxy);
                end

                % Smooth if needed
                if p.smooth
                    movslice = movmean(movslice, p.smoothwin * p.fps, 3);
                end

                % Movmed if needed
                if p.movmed
                    movslice = movmedian(movslice, p.movmedwin * p.fps, 3);
                end

                % Remove trim
                if p.smooth || p.movmed
                    movslice = movslice(:,:,trim+1 : end-trim);
                end

                % Bint if needed
                if p.bint > 1
                    movslice = bint(movslice, p.bint);
                end
                
                if p.dffbeforeaverage && p.dffmovie
                    % Average frame
                    movslice_base = mean(movslice(:,:,1:dffframes2), 3);
                    movslice = movslice ./ repmat(movslice_base, [1 1 l2]) - 1;
                end
                
                if p.znorm
                    movslice = movslice / zmf(k) * zmf_m;
                end

                % Save
                mov4d_f_cell{cellid_curr, ifov}(:,:,:,fi) = movslice;
            end
        end
    end
end

%% Post processing
% Keep track of cell count
ind_total_l = 0;
ind_total_r = 0;

% Canvas size
canvassize = panelsize2.* grid;

% Initialize canvas
canvasl = zeros(canvassize(1),canvassize(2),l2);
canvasr = zeros(canvassize(1),canvassize(2),l2);

% Keep track of the lowest row
gridy_max = 0;

if ~p.dffmovie
    canvasl = uint16(canvasl);
    canvasr = uint16(canvasr);
end

for i = 1 : nfov
    for j = 1 : mcells
        % Current np and nf
        np = np_mat(j, i);
        nf = nf_mat(j, i);
                
        % Do the averaging between trials
        if np > 0
            mov4d_p = squeeze(mean(mov4d_p_cell{j, i},4));
            ind_total_l = ind_total_l + 1;
        end
        if nf > 0
            mov4d_f = squeeze(mean(mov4d_f_cell{j, i},4));
            ind_total_r = ind_total_r + 1;
        end
        
        if strcmp(p.tiling, 'matching')
            ind_total_l = max(ind_total_l, ind_total_r);
            ind_total_r = max(ind_total_l, ind_total_r);
        end
        
        % dff
        if p.dffmovie
            if ~p.dffbeforeaverage
                % Average frame
                if np > 0
                    mov4d_p_base = mean(mov4d_p(:,:,1:dffframes2), 3);
                end
                if nf > 0
                    mov4d_f_base = mean(mov4d_f(:,:,1:dffframes2), 3);
                end

                % dff movie
                if np > 0
                    mov4d_p = mov4d_p ./ repmat(mov4d_p_base, [1 1 l2]) - 1;
                end
                if nf > 0
                    mov4d_f = mov4d_f ./ repmat(mov4d_f_base, [1 1 l2]) - 1;
                end
            end

            % White box
            if p.usewhitebox
                if np > 0
                    mov4d_p(1:p.whiteboxsize, 1:p.whiteboxsize, preframes2+1 : preframes2+boxlength)...
                        = p.whiteboxvalue(2);
                end
                if nf > 0
                    mov4d_f(1:p.whiteboxsize, 1:p.whiteboxsize, preframes2+1 : preframes2+boxlength)...
                        = p.whiteboxvalue(2);
                end
            end
            
            % Mask
            if p.showmask
                if np > 0 || nf > 0
                    % Make mask
                    mask = masks_xrun{j, i};
                    mask = repmat(mask, [1 1 l2]);
                    
                    % Label mask
                    if np > 0
                        mov4d_p(mask) = p.maskvalue(2);
                    end
                    if nf > 0
                        mov4d_f(mask) = p.maskvalue(2);
                    end
                end
            end
        else
            % White box
            if p.usewhitebox
                if np > 0
                    mov4d_p(1:p.whiteboxsize, 1:p.whiteboxsize, preframes2+1 : preframes2+boxlength)...
                        = p.whiteboxvalue(1);
                end
                if nf > 0
                    mov4d_f(1:p.whiteboxsize, 1:p.whiteboxsize, preframes2+1 : preframes2+boxlength)...
                        = p.whiteboxvalue(1);
                end
            end
            
            % Mask
            if p.showmask
                if np > 0 || nf > 0
                    % Make mask
                    mask = masks_xrun{j, i};
                    mask = repmat(mask, [1 1 l2]);
                    
                    % Label mask
                    if np > 0
                        mov4d_p(mask) = p.maskvalue(1);
                    end
                    if nf > 0
                        mov4d_f(mask) = p.maskvalue(1);
                    end
                end
            end
            
            % Turn into integer
            if np > 0
                mov4d_p = uint16(mov4d_p);
            end
            if nf > 0
                mov4d_f = uint16(mov4d_f);
            end
        end
        
        % Save to canvas
        if np > 0
            % Find location on the grid
            gridx = mod(ind_total_l, p.nfperrow);
            if gridx == 0
                gridx = p.nfperrow;
            end

            if gridx == p.nfperrow
                gridy = floor(ind_total_l/p.nfperrow);
            else
                gridy = floor(ind_total_l/p.nfperrow) + 1;
            end
            canvasl((gridy-1) * panelsize2(1)+1 : gridy * panelsize2(1),...
                (gridx-1) * panelsize2(2)+1 : gridx * panelsize2(2), :) = mov4d_p;
            
            gridy_max = max(gridy_max, gridy);
        end
        
        
        if nf > 0
            % Find location on the grid
            gridx = mod(ind_total_r, p.nfperrow);
            if gridx == 0
                gridx = p.nfperrow;
            end

            if gridx == p.nfperrow
                gridy = floor(ind_total_r/p.nfperrow);
            else
                gridy = floor(ind_total_r/p.nfperrow) + 1;
            end
            canvasr((gridy-1) * panelsize2(1)+1 : gridy * panelsize2(1),...
                (gridx-1) * panelsize2(2)+1 : gridx * panelsize2(2), :) = mov4d_f;
            gridy_max = max(gridy_max, gridy);
        end
    end
end
    
% Concatenate
if p.separationbar
    if p.dffmovie
        barvalue = p.whiteboxvalue(2);
    else
        barvalue = p.whiteboxvalue(1);
    end
    
    if p.concatenatedir == 2
        bar = ones(canvassize(1), p.separationbarsize, l2) * barvalue;
    else
        bar = ones(p.separationbarsize, canvassize(2), l2) * barvalue;
    end
    canvas = cat(p.concatenatedir, canvasl, bar, canvasr);
else
    canvas = cat(p.concatenatedir, canvasl, canvasr);
end


% trim
if gridy_max < grid(1)
    % Trim
    canvas = canvas(1 : panelsize2(1) * gridy_max, :,:);
end
%% Output
% Movie
writetiff(canvas, fullfile(p.outputpath, p.outputfn));

% Paramters
if p.savepars
    pars = struct;
    pars.mice = mousecell;
    pars.dates = datecell;
    pars.runs = runvec;
    pars.OTs = optotunecell;
    pars.p = p;
    pars.np_mat = np_mat;
    pars.nf_mat = nf_mat;
    pars.fovvec = fovvec;
    pars.ip_mat = ip_mat;
    pars.if_mat = if_mat;
    pars.passindcell = passindcell;
    pars.failindcell = failindcell;
    pars.cellid_xrun = cellid_xrun;
    pars.ncell_xrun = ncell_xrun;
    pars.xrunstrs = xrunstrs;

    outputfn_par = sprintf('%s.mat', p.outputfn(1:end-4));
    save(fullfile(p.outputpath, outputfn_par), '-struct', 'pars', '-v7.3');
end

% Materials
if p.savematerials
    outputfn_material = sprintf('%s_materials.mat', p.outputfn(1:end-4));
    save(fullfile(p.outputpath, outputfn_material), 'mov4d_p_cell', 'mov4d_f_cell', '-v7.3');
end


end