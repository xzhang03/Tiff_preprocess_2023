function TiffCombineOptoSplits(varargin)

if nargin < 1
    varargin = {};
end

p = inputParser;
addOptional(p, 'defaultfp', '\\nasquatch\data\2p\stephen\'); 
addOptional(p, 'defaultfn', '*.tif');

parse(p, varargin{:});
p = p.Results;

%% IO
[fns, fp, ~] = uigetfile(fullfile(p.defaultfp, p.defaultfn), 'Select the tiff files.', ...
    'MultiSelect', 'on');
n = length(fns);

%% Figure out how to combine
switch n
    case 2
        r = 1;
        c = 2;
    case 4
        r = 2;
        c = 2;
    case 8
        r = 2;
        c = 4;
    case 9
        r = 3;
        c = 3;
    case 10
        r = 2;
        c = 5;
    case 15
        r = 3;
        c = 5;
    case 20
        r = 4;
        c = 5;
    case 25
        r = 5;
        c = 5;
    case 30
        r = 5;
        c = 6;
    case 35
        r = 5;
        c = 7;
    case 40
        r = 5;
        c = 8;
    otherwise
        r = [];
        c = [];
        
end
definput = {num2str(r), num2str(c), '1'};
answer = inputdlg({'How many rows:', 'How many columns:', 'Bin XY:'}, ...
    sprintf('How to combine %i movies', n),...
    [1 30], definput);
r = str2double(answer{1});
c = str2double(answer{2});
bin_xy = str2double(answer{3});

%% Load
movs = cell(n,1);
hwait = waitbar(0, 'Processing');
for i = 1 : n
    waitbar(i/n, hwait, sprintf('Loading %i/%i fies', i, n));
    movs{i} = readtiff(fullfile(fp, fns{i}));
    if bin_xy > 1
        movs{i} = binxy(movs{i}, bin_xy);
    end
end
close(hwait)

%% Combining horizontal
rs  = cell(r, 1);
ind = 0;
hwait = waitbar(0, 'Processing');
for i = 1 : r
    cmov = [];
    for j = 1 : c
        ind = ind + 1;
        waitbar(ind/n, hwait, sprintf('Horizontal combining %i/%i fies', ind, n));
        cmov = cat(2, cmov, movs{ind});
    end
    rs{i} = cmov;
end
close(hwait)

%% Combining vertical
outputmov = [];
hwait = waitbar(0, 'Processing');
for i = 1 : r
    waitbar(i/r, hwait, sprintf('Vertical combining %i/%i fies', i, r));
    outputmov = cat(1, outputmov, rs{i});
end
close(hwait)

%% Output
fnout = sprintf('%s%i-%i.tif', fns{1}(1:end-5), 1, n);
writetiff(outputmov, fullfile(fp, fnout), 'double');

end