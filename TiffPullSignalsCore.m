function cellsort = TiffPullSignalsCore(mov, cellsort, xybin)
%TiffPullSignalsCore pulls signals from tiff data

   
% Get original file size and number of frames
nrois = length(cellsort);
nframes = size(mov,3);

% Initialize
signal = zeros(size(mov, 3), nrois);
neuropil = zeros(size(mov, 3), nrois);

mov = reshape(mov, size(mov,1)*size(mov,2), nframes);

for j = 1:nrois
    % Get mask
    mask = cellsort(j).mask;
    mask = binxy(mask, xybin) > 0;
    
    % Get neuropil
    np = cellsort(j).neuropil;
    np = binxy(np, xybin) > 0;
    
    % Pull
    signal(:, j) = mean(mov(mask, :));
    neuropil(:, j) = mean(mov(np, :));
end

% Subtract
signalsub = signal - neuropil;


% Now put into cellsort format
for r = 1 : nrois
    cellsort(r).timecourse.raw = signal(:, r)';
    cellsort(r).timecourse.neuropil = neuropil(:, r)';
    cellsort(r).timecourse.subtracted = signalsub(:,r)';
end

end

