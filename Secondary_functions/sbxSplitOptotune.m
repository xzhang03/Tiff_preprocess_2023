function sbxSplitOptotune(path, k, N, pmt)
% sbxSplitOptotune Splits an sbx file into tiff files of the same optotune
% level
% sbxSplitOptotune(path, k, N, pmt)
%
% Reads from frame k to k + (N - 1) in file fname
% 
% path  - the file path to .sbx file (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.,
% optional
% pmt   - the number of the pmt, 0 for green or 1 for red, assumed to be 0

% If N>1 it returns a 4D array of size = [#pmt rows cols N] 
% If N=1 it returns a 3D array of size = [#pmt rows cols]
% If N<0 it returns an array to the end

% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
%
%
% The function also creates a global 'info' variable with additional
% informationi about the file
% 2021/05/18 Modified to work with optotune in scanbox-tower -SZ

% Force a reload of the global info variables. Without this, troube arises
%clearvars -global info 
info = sbxInfo(path, true);
% Check if optotune was used, accounting for the version of scanbox being
% used
optotune_used = false;
if isfield(info, 'volscan') && info.volscan > 0, optotune_used = true; end
if ~isfield(info, 'volscan')
    if isfield(info, 'otwave') && ~isempty(info.otwave)
        optotune_used = true;
    elseif isfield(info, 'etl_table') && size(info.etl_table,1) > 1
        optotune_used = true;
    end
end

% Set to start at beginning if necessary
if nargin < 2, k = 0; end
% Set in to read the whole file if unset
if nargin < 3 || N < 0, N = info.max_idx + 1 - k; end
% Make sure that we don't search beyond the end of the file
if N > info.max_idx + 1 - k, N = info.max_idx + 1 - k; end

% Automatically set the PMT to be 0
if nargin < 4, pmt = 0; end

% Fix 0-to-1 indexing
pmt = pmt + 1;

if ~optotune_used
    warning('No optotune detected');
    return
end

% Number of frames
if isfield(info, 'otwave')
    nOTlevels = length(info.otwave);
elseif isfield(info, 'etl_table')
    nOTlevels = size(info.etl_table,1);
end

fprintf('Found %i optotune levels.\n', nOTlevels);

% File names
[fp, fn_sbx, ~] = fileparts(path);
fns = cell(nOTlevels, 1);
for i = 1 : nOTlevels
    fn = fullfile(fp, sprintf('%s_OT%i-%i.tif', fn_sbx, i, pmt));
    if exist(fn, 'file')
        overwrite = input(sprintf('%s_OT%i.tif already exists! Overwrite (1 = yes, 0 = no)? ', fn_sbx, i));
        if overwrite ~= 1
            disp('Overwrite canceled.')
            return;
        end
    end
    fns{i} = fn;
end

if (isfield(info, 'fid') && info.fid ~= -1)
    try
        fseek(info.fid, k*info.nsamples, 'bof');
        x = fread(info.fid, info.nsamples/2*N, 'uint16=>uint16');
        x = reshape(x, [info.nchan info.sz(2) info.recordsPerBuffer N]);
    catch
        error('Cannot read frame. Index range likely outside of bounds.');
    end

    x = intmax('uint16') - permute(x, [1 3 2 4]);
    
    % Added by Arthur-- correct the output to a single PMT
    if info.nchan == 1
        if N > 1
            x = squeeze(x(1, :, :, :));
        else
            x = squeeze(x(1, :, :)); 
        end
    else
        if N > 1
            x = squeeze(x(pmt, :, :, :));
        else
            x = squeeze(x(pmt, :, :));
        end
    end
    
    % Splitting
    for i = 1 : nOTlevels
        % Split
        optoframes = i : nOTlevels : size(x,3);
        xot = x(:,:,optoframes);
        
        % Write
        writetiff(xot, fns{i});
    end
end