function TiffStrackReg(fpath,varargin)
%Registers a stack from center to the top/bottom

if nargin < 2
    varargin = {};
    if nargin < 1
        fpath = '';
    end
end

%% Parse inputs
p = inputParser;

% General path
addOptional(p, 'genpath', '\\nasquatch\data\2p\stephen');
addOptional(p, 'ext', '*.tif');

% Bin
addOptional(p, 'binxy', 2);

% Local Normalize
addOptional(p, 'localnormalize', true);
addOptional(p, 'lnvec', [8, 60]);

% Reg methods
addOptional(p, 'xyreg', true);
addOptional(p, 'demonsreg', true);

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;

if isempty(fpath)
    [fn, fp] = uigetfile(fullfile(p.genpath, p.ext));
    fpath = fullfile(fp, fn);
end

%% IO
% Load
mov = readtiff(fpath);
sz = size(mov);
n = size(mov, 3);

% Bin
movb = binxy(mov, p.binxy);
if p.localnormalize
    for i = 1 : n
        movb(:,:,i) = localnormalizecore(movb(:,:,i), p.lnvec);
    end
end

%% Reg
% Start frame
istart = round(n/2);

% Intialize
movreg = mov;
movreg(:,:,istart) = mov(:,:,istart);

% Get the meshgrid to pass to all the workers
[xx, yy] = meshgrid(1 : sz(2), 1 : sz(1));

% Reg to the top
hwait = waitbar(0);
f = 0;
for i = istart-1 : -1 : 1
    % Progress
    f = f + 1;
    waitbar(f/n, hwait, sprintf('Registering to front %i/%i', f, n));
    
    % XY reg
    if p.xyreg
        [shifts,~]=stackRegisterMA_RR(movb(:,:,i), movb(:,:,i+1), 100, [], 0);
        [~, movb(:,:,i)] = stackRegisterMA_RR(movb(:,:,i), [], [], shifts, 0);
        [~, movreg(:,:,i)] = stackRegisterMA_RR(mov(:,:,i), [], [], shifts * p.binxy, 0);
    end
    
    % Demons reg
    if p.demonsreg
        [D,~] = imregdemons(movb(:,:,i), movb(:,:,i+1), [32 16 8 4],...
                'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitbar',false);
        D = imresize(D, p.binxy);
        movreg(:,:,i) = interp2(xx, yy, movreg(:,:,i), xx + D(:,:,1), yy + D(:,:,2));
    end
end

% Reg to the bottom
for i = istart+1 : 1 : n
    % Progress
    f = f + 1;
    waitbar(f/n, hwait, sprintf('Registering to back %i/%i', f, n));
    
    % XY reg
    if p.xyreg
        [shifts,~]=stackRegisterMA_RR(movb(:,:,i), movb(:,:,i-1), 100, [], 0);
        [~, movb(:,:,i)] = stackRegisterMA_RR(movb(:,:,i), [], [], shifts, 0);
        [~, movreg(:,:,i)] = stackRegisterMA_RR(mov(:,:,i), [], [], shifts * p.binxy, 0);
    end
    
    % Demons reg
    if p.demonsreg
        [D,~] = imregdemons(movb(:,:,i), movb(:,:,i-1), [32 16 8 4],...
                'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitbar',false);
        D = imresize(D, p.binxy);
        movreg(:,:,i) = interp2(xx, yy, movreg(:,:,i), xx + D(:,:,1), yy + D(:,:,2));
    end
end
close(hwait)

%% Remove nan
for i = 1 : n
    f = movreg(:,:,i);
    if any(isnan(f(:)))
        nanv = nanmean(f(:));
        f(isnan(f)) = nanv;
        movreg(:,:,i) = f;
    end
end

%% Write
[fp, fn, ext] = fileparts(fpath);
fnout = sprintf('%s_streg%s', fn, ext);
writetiff(movreg, fullfile(fp, fnout), 'double');

end

