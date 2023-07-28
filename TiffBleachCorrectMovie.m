function movout = TiffBleachCorrectMovie(mov, varargin)
%TiffBleachCorrectMovie correct bleaching on a pixel level.

% Parse inputs
p = inputParser;

% Function
addOptional(p, 'fun', 'a*exp(-b*x)+c');
addOptional(p, 'lowerbound', [-Inf, 0, 0]);
addOptional(p, 'startingpoint', [1000 0.0001 1000]);

% Data points
addOptional(p, 'seg1', [1 1000]);
addOptional(p, 'seg2', []); % Leave empty to correct only based on segment 1. Careful.

% Output
addOptional(p, 'savetiff', false);
addOptional(p, 'savepath', '');

% Unpack if needed
if size(varargin,1) == 1 && size(varargin,2) == 1
    varargin = varargin{:};
end

% Parse
parse(p, varargin{:});
p = p.Results;

type = class(mov);

%% Get the frames
% fitting frames
fs1 = (p.seg1(1) : p.seg1(2))';
if ~isempty(p.seg2)
    fs2 = (p.seg2(1) : p.seg2(2))';
else
    fs2 = [];
end
fs = cat(1, fs1, fs2);

% all frames
sizevec = size(mov);
l = sizevec(3);

%% Fit overall
% Overall
overall = squeeze(mean(mean(mov,1),2));

% Exp fitting
f1_exp = fit(fs, overall(fs), 'a*exp(-b*x)+c',...
    'Lower', [-Inf, 0, 0], 'StartPoint', [0.2, 0.0001, 4]);

% Matrices for linear fitting
f1 = f1_exp(fs);
f2 = ones(length(f1),1);
f12 = [f1 f2];

f1_long = f1_exp(1:l);
f2_long = ones(length(f1_long),1);
f12_long = [f1_long f2_long];

%% Fit pixel by pixel
% Initialize
movout = zeros(size(mov));
movout = cast(movout, type);

tic
hwait = waitbar(0, 'Bleach correcting...');
for i = 1 : sizevec(1)
    if mod(i,10) == 0
        waitbar(i/sizevec(1), hwait, sprintf('Bleach correcting %i/%i', i, sizevec(1)));
    end
    
    for j = 1 : sizevec(2)
        vec = double(squeeze(mov(i, j,:)));
        s = f12 \ vec(fs);
        vec2 = vec - f12_long * s;
        vec2 = vec2 + mean(vec(fs1)) - mean(vec2(fs1));
        movout(i,j,:) = vec2;
    end
end
close(hwait)
t = round(toc);
fprintf('Bleaching correction done. Elapsed time = %i seconds.\n', t);

%% Save
if p.savetiff
    writetiff(movout, p.savepath);
end
end