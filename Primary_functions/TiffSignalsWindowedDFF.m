function cellsort = TiffSignalsWindowedDFF(cellsort, fps, time_window, percentile)
% TiffSignalsWindowedDFF Get the dff trace and axon dff trace from a cellsort and
%   add it to the cellsort. This is the non-parallel version of the code

    % time window: moving window of X seconds - calculate f0 at time window prior to each frame - used to be 32
    if nargin < 3, time_window = 32; end
    if nargin < 4, percentile = 10; end

    nframes = length(cellsort(1).timecourse.raw);
    nrois = length(cellsort);

    % Calculate f0 for each timecourse using a moving window of time window
    % prior to each frame
    f0 = zeros(nrois, nframes);
    windowframes = round(time_window*fps);

    % create temporary traces variable that allows you to do the prctile
    % quickly
    traces_f = nan(nrois, length(cellsort(1).timecourse.subtracted));
    hwait = waitbar(0, sprintf('DFF ROI: %i/%i', 0, nrois));
    for curr_ROI = 1:nrois
        if mod(curr_ROI, 10) == 0
            waitbar(curr_ROI/nrois, hwait, sprintf('DFF ROI: %i/%i', curr_ROI, nrois))
        end
        traces_f(curr_ROI, :) = cellsort(curr_ROI).timecourse.subtracted;
        
        for i = 1:nframes
            if i <= windowframes
                frames = traces_f(curr_ROI,1:windowframes);
                f0(curr_ROI,i) = prctile(frames,percentile);
            else
                frames = traces_f(curr_ROI,i - windowframes:i-1);
                f0(curr_ROI,i) = prctile(frames,percentile);
            end
        end
    end
    close(hwait)

    % dff
    traces_dff = (traces_f-f0)./ f0;

    % Stick back into cellsort variable
    for curr_ROI = 1:nrois
        cellsort(curr_ROI).timecourse.f0_axon = f0(curr_ROI,:);
        cellsort(curr_ROI).timecourse.dff_axon = traces_dff(curr_ROI,:);
        cellsort(curr_ROI).timecourse.dff_axon_norm = cellsort(curr_ROI).timecourse.dff_axon ./ max(cellsort(curr_ROI).timecourse.dff_axon);
    end    
end

