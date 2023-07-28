%% Identifiers and edges
clear
mouse = 'SZ871';
date = '220216';
optoset = [1, 2];
runs = [1];
fps = 30.77;

% Edge
edges = [124,100,24,8];

%% Split
sbxSplitOptotune(datapath, 0, -1, 0);

%% Standard xyreg
for i = 1 : optoset(2)
    Tiffxyreg(mouse, date, runs, 'server', 'nasquatch', 'pmt', 1, 'optotune', i, 'edges', edges,...
        'chunksize', 800, 'binxy', 2, 'gausssize', [16 60], 'reuseref', false, 'autofocus', true);
end

%% Demonsereg
for i = 1 : optoset(2)
    sbxDemonsRegOneshotTiff(mouse, date, 'runs', runs, 'server', 'nasquatch', 'pmt', 1,...
        'optotune', i, 'edges', edges, 'chunksize', 100, 'binxy', 2, 'hp_norm_sigmas', [16 60],...
        'binbeforehighpassnorm', true, 'movtype', 'OTtiff_xyreg', 'savewarp', false, 'reuseref', true,...
        'AccumulatedFieldSmoothing', 1, 'parfor', true, 'nworkers', 4);
end

%% Morphological filter
% Just used to generate segmentation file now
tic
for i = 1 : optoset(2)
    TiffMorphologicalFilterExtractROIs(mouse, date, runs(1), 'server', 'nasquatch', 'optotune', i,...
        'downsample_xy', 2, 'movtype', 'OTtiff_demonsreg', 'threshold', 0.05, 'justmeanim', true);
end
toc

%% Cellpose then
% Use cellpose
return;

%% Signals
for i = 1 : optoset(2)
    TiffSignalsCellpose(mouse, date, runs, 'optotune', i, 'freq', 10, 'server', 'nasquatch', 'pmt', 1,...
        'optotune', i, 'binxy', 2, 'freq', fps/optoset(2), 'percentile', 90, 'edges', [], 'GRIN', true);
end

