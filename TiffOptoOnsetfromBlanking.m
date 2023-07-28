function [nons, noffs, vector] = TiffOptoOnsetfromBlanking(im, varargin)
% Get opto onset from pmt blanking (top of hte screen)

if nargin < 2
    varargin = {};
end

% Parse
% Whether or not axons was used is determined from icaguidata.pars
p = inputParser;
addOptional(p, 'rows', [1 10]);  % PMT to use for extraction

parse(p, varargin{:});
p = p.Results;

%% IO
if ischar(im)
    im = readtiff(im);
end

%% Get trace
% Calculate trace
t = squeeze(nanmean(nanmean(im(p.rows(1) : p.rows(2), :, :),1),2));

%% Draw bbox
% Draw
hfig = figure;
plot(t)
title('Draw a bounding box for blanked time and intensity range');
r = imrect(); 
bbox = wait(r);
close(hfig)

% Bounds
bounds = round([bbox(1), bbox(1)+bbox(3), bbox(2), bbox(2)+bbox(4)]);

%% Determine points
% Vector form
vector = (t >= bounds(3)) & (t <= bounds(4));
vector(1 : bounds(1)-1) = false;
vector(bounds(2)+1 : end) = false;

% Mat form
chains = chainfinder(vector);

% Output
nons = chains(:,1);
noffs = chains(:,1) + chains(:,2) - 1;

end

