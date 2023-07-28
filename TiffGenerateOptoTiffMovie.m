function movout = TiffGenerateOptoTiffMovie(mouse,date,run, pmt,sbxtype,varargin)
%Expects trialtype to be a cell array of cell arrays
p = inputParser;
addOptional(p, 'fr', []); %allows for manual entering frame rate
addOptional(p, 'prestim', 2); % in seconds
addOptional(p, 'poststim', 8); % in seconds. Can leave 0 to make the video continuous.
addOptional(p, 'bin_xy', 2);
addOptional(p, 'bin_t', 2);
addOptional(p, 'given_nons', []); %allows manually providing frame onsets of interest
addOptional(p, 'given_noffs', []); %allows manually providing frame offsets of interest
addOptional(p, 'noduplicate', false);
addOptional(p, 'save', true);
addOptional(p, 'passmovie', false);
addOptional(p, 'server', []);
addOptional(p, 'framechannel', 2); % Frame channel in the ephys file
addOptional(p, 'optochannel', 8); % Where to find opto TTL on nidaq
addOptional(p, 'optotraingap', 5); % Gaps of trains of opto stims. In seconds.
addOptional(p, 'trialtype', ''); % A string for naming purpose only
addOptional(p, 'XZmovie', []);

% Stim indicator
addOptional(p, 'usewhitesquare', true);
addOptional(p, 'multiwhitesquare', true);

% Movie type
addOptional(p, 'fromraw', false); % Raw or dff movie
addOptional(p, 'makeaveragestimmovie', false);
addOptional(p, 'makemoviesplits', false);

% External trigger
addOptional(p, 'useexternaltrig', false); % Use external trigger
addOptional(p, 'min_frac_amplitude', 0.5); % Minimal fraction drop amplitude
addOptional(p, 'polarity', 'neg'); % Polarity of the external trigger
addOptional(p, 'forceonduration', 0); % Number of frames to force on the on-trigger, set to 0 to turn off

% Premedian filter
addOptional(p, 'premedianfilter', []);

% Prebeach correct
addOptional(p, 'prebleachcorrect', false);
addOptional(p, 'fun', 'a*exp(-b*x)+c');
addOptional(p, 'lowerbound', [-Inf, 0, 0]);
addOptional(p, 'startingpoint', [1000 0.0001 1000]);
addOptional(p, 'seg1', [1 1000]);
addOptional(p, 'seg2', []); % Leave empty to correct only based on segment 1. Careful.

% Trials
addOptional(p, 'trials', []); % Leave blank to use all trials

% Optotune
addOptional(p, 'optotune', []); % Optotune

% Debug
addOptional(p, 'mov', []);% Directly pass movie here
addOptional(p, 'rev', false);

parse(p, varargin{:});
p = p.Results;

% Set time in seconds before and after stim
prestim = p.prestim;
poststim = p.poststim;

% % Get the onset times
% ons = sbxOnsets(mouse, date, run, p.server);
% if isfield(ons, 'onsets') < 1, return; end

% Load in the movies
infopath = sbxPath(mouse, date, run, 'sbx', 'pmt',pmt,'server',p.server,'estimate',false); % SZ changed estimate to false 4/3/2019
info = sbxInfo(infopath);
fr = 30.98;
if info.scanmode > 0, fr = 15.49; end

if ~isempty(p.fr); fr = p.fr; end

if ~isempty(p.optotune)
    noptotune = size(info.etl_table,1);
    fr = fr / noptotune;
end

if ~isempty(p.given_nons)
    nons = p.given_nons; 
    noffs = p.given_noffs;
elseif p.useexternaltrig
    % Read
    [fp, ~, ~] = fileparts(movpath);
    [fn, fp]= uigetfile(fullfile(fp, '*.csv'));
    exttrigger = readtable(fullfile(fp,fn));
    exttrigger = table2array(exttrigger(:,end));
    
    % Threshold
    thresh = median(exttrigger) * p.min_frac_amplitude;
    
    switch p.polarity
        case 'pos'
            nons = find(diff(exttrigger) >= thresh) + 1;
            if p.forceonduration == 0
                noffs = find(diff(exttrigger) <= -thresh) + 1;
            else
                noffs = nons + p.forceonduration;
            end
        case 'neg'
            nons = find(diff(exttrigger) <= -thresh) + 1;
            if p.forceonduration == 0
                noffs = find(diff(exttrigger) >= thresh) + 1;
            else
                noffs = nons + p.forceonduration;
            end
    end
else
    [nons, noffs] = sbxOnsetsfromNidaq(mouse, date, run, 'server', p.server, ...
        'pulsechannel', p.optochannel, 'traingap', p.optotraingap, 'framechannel',...
        p.framechannel);
end

if ~isempty(p.trials)
    nons = nons(p.trials);
    noffs = noffs(p.trials);
end

% Adjust for optotune
if ~isempty(p.optotune)
    nons = round(nons / noptotune);
    noffs = round(noffs / noptotune);
end

% Frames before stim
preframes = ceil(prestim*fr);

% Frames after stim
if isempty(poststim)
    % Continuous mode
    % Get a concensus ITI
    ITI = mode(diff(nons));
    postframes = ITI - preframes; 
else
    postframes = ceil(poststim*fr);
end

% Movie path
if isempty(p.optotune)
    movpath = sbxPath(mouse, date, run, sbxtype, 'pmt',pmt,'server',p.server,'estimate',false); % SZ changed estimate to false 4/3/2019
else
    movpath = sbxPath(mouse, date, run, sbxtype, 'pmt',pmt,'server',p.server,'estimate',false, 'optotune', p.optotune); % SZ changed estimate to false 4/3/2019
end

if ~isempty(nons)
    if p.fromraw
        mt = 'raw';
    else
        mt = 'dff';
    end
    
    spath = sprintf(['%s_%s-%s-across-%i_from' sbxtype '.tif'], movpath(1:strfind(movpath, '.')-1), mt, p.trialtype, length(nons));

    nim =[];
    trialmov =[];
    
    % White square info
    if p.usewhitesquare
        wsframes = zeros(length(nons), 2);
        for i = 1 : length(nons)
            wsframes(i, 1) = preframes + 1;
            wsframes(i, 2) = preframes + noffs(i) - nons(i) + 1;
        end
        
        % Multi whitesquares
        if p.multiwhitesquare
            dnons = diff(nons);
            dnons = dnons + noffs(2:end) - nons(2:end);
            nws = floor(postframes ./ dnons) + 1;
            nws(end+1) = 1;
            
            for i = 1 : length(nons)
                for ii = 2 : nws(i)
                    if i+ii-1 <= length(nons)
                        wsframes(i, (ii-1)*2+1) = wsframes(i, (ii-2)*2+1) + nons(i+ii-1) - nons(i+ii-2);
                        wsframes(i, (ii-1)*2+2) = wsframes(i, (ii-2)*2+2) + noffs(i+ii-1) - noffs(i+ii-2);
                    end
                end
            end
        else
            nws = ones(length(nons),1);
        end
    end
    
    % Remove duplicate
    if p.noduplicate
        nwsmax = max(nws);
        nons = nons(1:nwsmax:end);
        noffs = noffs(1:nwsmax:end);
        wsframes = wsframes(1:nwsmax:end,:);
    end
    
    % Read movie
    if isempty(p.mov)
        mov = readtiff(movpath);
    else
        mov = p.mov;
        p = rmfield(p, 'mov');
    end
    
    % Bin once for all    
    if p.bin_xy > 1
        tic
        mov = binxy(mov, p.bin_xy);
        t = round(toc);
        fprintf('Binning xy done. Elapsed time = %i seconds.\n', t);
    end
    
    
    
    if p.bin_t > 1
        tic
        mov = bint(mov, p.bin_t);
        preframes = floor(preframes / p.bin_t);
        postframes = floor(postframes / p.bin_t);
        
        nons = round(nons / p.bin_t);
        noffs = round(noffs / p.bin_t);
        
        wsframes = round(wsframes / p.bin_t);
        
        p.seg1 = round(p.seg1 / p.bin_t);
        p.seg2 = round(p.seg2 / p.bin_t);
        t = round(toc);
        fprintf('Binning t done. Elapsed time = %i seconds.\n', t);
    end
    
    
    if ~isempty(p.premedianfilter)
        tic
        mov = movmedian(mov, p.premedianfilter, 3);
        t = round(toc);
        fprintf('Movmedian done. Elapsed time = %i seconds.\n', t);
    end
    
    if p.prebleachcorrect
        mov = TiffBleachCorrectMovie(mov, 'seg1', p.seg1, 'seg2', p.seg2,...
            'fun', p.fun, 'lowerbound', p.lowerbound, 'startingpoint', p.startingpoint);
    end
    
    tic;
    mov = double(mov);
    t = round(toc);
    fprintf('Conversion done. Elapsed time = %i seconds.\n', t);
        
    hwait = waitbar(0, 'Processing');
    for j = 1:length(nons)
        waitbar(j/length(nons), hwait, sprintf('Processing sweep %i/%i.', j, length(nons)))
        
        if nons(j) + postframes > size(mov,3)
            fprintf('Trial clipped. Ending at onset: %i\n', j);
            break
        end
        if p.fromraw
            trialmov = mov(:,:,nons(j) - preframes : nons(j) + postframes);
        else
            denom = averaget(mov(:,:,nons(j) - preframes : nons(j)));
            
            trialmov = mov(:,:,nons(j) - preframes : nons(j) + postframes);
            
            denom2 = repmat(denom,1,1,size(trialmov,3));
            
            trialmov = (trialmov - denom2) ./ denom2;
        end
        
        if p.rev
            trialmov = trialmov(:,:,end:-1:1);
        end
        
        if p.usewhitesquare
            if j == 1 % First time of stim, make white square
                redowhitesquare = true;
            elseif (noffs(j) - nons(j)) ~= (noffs(j-1) - noffs(j-1)) % Stim size changed, redo white square
                redowhitesquare = true;
            else
                redowhitesquare = false;
            end
        end
        
        if redowhitesquare
            whitesquaremov = zeros(size(trialmov,1),size(trialmov,2),preframes+postframes+1);
            if p.fromraw
                ws_value = 65534;
            else
                ws_value = 20;
            end
            
            for ii = 1 : nws(j)
                wsstart = wsframes(j, (ii-1)*2+1);
                wsend = wsframes(j, ii*2);
                
                if wsstart > 0 && wsend > 0
                    whitesquaremov(1:round(.1*size(trialmov,1)),1:round(.1*size(trialmov,2)),...
                        wsstart : wsend) = ws_value;
                end
            end
        end
        
        trialmov(isnan(trialmov)) = 0;
        if p.usewhitesquare
            trialmov = trialmov + whitesquaremov;
        end
        
        if ~isempty(nim)
            if p.makeaveragestimmovie || p.makemoviesplits
                nim = cat(1,nim,reshape(trialmov,1,size(trialmov,1)*size(trialmov,2)*size(trialmov,3)));
            else
                nim = cat(3,nim,trialmov);
            end
        else
            if p.makeaveragestimmovie || p.makemoviesplits
                nim = reshape(trialmov,1,size(trialmov,1)*size(trialmov,2)*size(trialmov,3));
            else
                nim = trialmov;
            end
        end
    end
    close(hwait)
    
    % debug
%     nim2 = reshape(nim, [length(nons) 256 398 724]);
%     nim2_a = mean(nim2, 2);
%     nim2_a = mean(nim2_a, 3);
%     nim2_a = squeeze(nim2_a);
    
    if p.makeaveragestimmovie
        nim = reshape(mean(nim, 1),size(trialmov));
        spath = sprintf('%s_MeanTrialMovie.tif', spath(1:strfind(spath, '.')-1));
    elseif p.makemoviesplits
        [fp, fn, ~] = fileparts(movpath);
        foldername = sprintf('%s_stimsplits', mouse);
        fp_split = fullfile(fp,foldername);
        if isdir(fp_split)
            mkdir(fp_split);
        end
    end
    
    if ~isempty(p.XZmovie)
        nim = nim(:,:,1:floor(size(nim,3)/size(info.otwave,2))*size(info.otwave,2));
        nimZ = reshape(nim,size(nim,1),size(nim,2),size(info.otwave,2),[]);
        
        if p.save
            writetiff(squeeze(nimZ(:,p.XZmovie,:,:)), spath, 'double');
        end
        if p.passmovie
            movout = squeeze(nimZ(:,p.XZmovie,:,:));
        else
            movout = [];
        end
    else
        if p.save
            if p.makemoviesplits
                hwait = waitbar(0, 'Processing');
                for j = 1 : length(nons)
                    waitbar(j/length(nons), hwait, sprintf('Saving sweep %i/%i.', j, length(nons)))
                    fnout = sprintf('%s_%s_Trial_%i.tif', fn, mt, j);
                    nim2 = reshape(nim(j,:),size(trialmov));
                    writetiff(nim2, fullfile(fp_split, fnout), 'double');
                end
                close(hwait)
            else
                writetiff(nim, spath, 'double');
            end
        end
        if p.passmovie
            movout = nim;
        else
            movout = [];
        end
    end
end
end