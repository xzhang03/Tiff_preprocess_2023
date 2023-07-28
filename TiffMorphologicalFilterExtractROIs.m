function TiffMorphologicalFilterExtractROIs(mouse, date, run, varargin)
%A tiff based morphological filter extraction code. Modified from previous
%andermann lab code
% TiffMorphologicalFilterExtractROIs(mouse, date, run, varargin)


p = inputParser;

% Lab variables variables
addOptional(p, 'server', []);  % Add in the server name as a string
addOptional(p, 'pmt', 1, @isnumeric);  % Which PMT to use for analysis, 1-green, 2-red
addOptional(p, 'filepath', ''); % Feed tif path directly

% Optotune level
addOptional(p, 'optotune', []); % Optotune level

% Frames for mean
addOptional(p, 'movsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
addOptional(p, 'movoffset', 500, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames
   
% IO variables
addOptional(p, 'mov', []); % Directly feed movie (debug)
addOptional(p, 'justmeanim', false); % Just make the mean im and leave (debug)
addOptional(p, 'movtype', 'OTtiff_demonsreg');  % input type, can be xyreg, sbxreg, or sbx.

addOptional(p, 'downsample_xy', 1);  % The half-width of cells in pixels
addOptional(p, 'threshold', 0.05);  %threshold for binarizing image
addOptional(p, 'n', 8);  %threshold for binarizing image
addOptional(p, 'm', 150);  %threshold for binarizing image
parse(p, varargin{:});
p = p.Results;

%% IO
% Path
if ~isempty(p.filepath)
    tiffpath = p.filepath;
else
    tiffpath = sbxPath(mouse, date, run, p.movtype, 'server', p.server, 'pmt', p.pmt, 'optotune', p.optotune);
end

% icapath
[outputfp, outputfn, ~] = fileparts(tiffpath);
if isempty(p.optotune)
    icapath = fullfile(outputfp, sprintf('%s_%s_%03d.ica', mouse, date, run));
else
    icapath = fullfile(outputfp, sprintf('%s_%s_%03d_OT%i.ica', mouse, date, run,p.optotune));
end

% meantifpath
if ~isempty(p.filepath)
    meanpath = fullfile(outputfp, sprintf('%s_toseg.tif', outputfn));
elseif isempty(p.optotune)
    meanpath = fullfile(outputfp, sprintf('%s_%s_%03d_toseg.tif', mouse, date, run));
else
    meanpath = fullfile(outputfp, sprintf('%s_%s_%03d_OT%i_toseg.tif', mouse, date, run,p.optotune));
end

% Read
if isempty(p.mov)
    mov = readtiff(tiffpath);
else
    mov = p.mov;
end

meanim = median(mov(:,:,p.movoffset : p.movoffset + p.movsize - 1),3);
meanim = double(meanim);
meanim = imresize(meanim,1/p.downsample_xy);

% Prepare target by local normalizing it
n=p.n; %assumes neurons are around 8 pixels wide
m=p.m;
    
meanim = medfilt2(double(meanim), [2, 2], 'symmetric');
meanim_prime = meanim - imgaussfilt(double(meanim),n);
ln_meanim = meanim_prime ./ (imgaussfilt(meanim_prime.^2,m) .^ (1/2));
ln_meanim(isnan(ln_meanim)) = 0;
    
writetiff(ln_meanim, meanpath);

if p.justmeanim
    return
end

subplot(1,2,1) 
imshow(meanim,[])
subplot(1,2,2) 
imshow(ln_meanim,[])

%%
Images = meanim;
pos = [];
Threshold = p.threshold;

Combined=zeros(size(Images(:,:,1)));
Combined2=zeros(size(Images(:,:,1)));
CountImage=zeros(size(Images(:,:,1)));

for(im=1:size(Images,3))
    f=Images(:,:,im);
    f(isnan(f))=0;
    f(isinf(f))=0;

    f=f-imopen(f,strel('rectangle',[18,3]));
    f=f-imopen(f,strel('rectangle',[3,18]));
    
    f_prime=double(f)-double(imgaussfilt(double(f),n));
    g=f_prime./(imgaussfilt(f_prime.^2,m).^(1/2));
    
    g(isnan(g))=0;
    g(isinf(g))=0;
    
    closed=imclose(g,strel('rectangle',[4,4]));
    opened=imopen(closed,strel('rectangle',[4,4]));
    
    %For Andrew (remove pipette)
    if ~isempty(pos)
        opened(pos(2):pos(4),pos(1):pos(3))=0;
    end
    
    % Remove borders
    opened(1:5,:)=0;
    opened(:,1:5)=0;
    opened((end-4:end),:)=0;
    opened(:,(end-4:end))=0;
    
    X=imdilate(opened,strel('disk',1));
    threshold_image=im2bw(double(X)/double(max(max(X))), Threshold);
    X=imdilate(threshold_image,strel('rectangle',[2,2]));
    A=bwdist(~X);
    B=A;
    B(A<=1.5)=0;
    Combined=Combined+B;
    A=-A;
    A(~X) = -Inf;
    A=watershed(A)-1;
    A=ceil(double(A)/max(max(double(A))));
    pos2 = pos;
    if ~isempty(pos);
    if im == 1;
        if pos(1)>5; pos2(1) = pos(1)-5; end
        if pos(2)>5; pos2(2) = pos(2)-5; end
        if pos(3)>5; pos2(3) = pos(3)+5; end
        if pos(4)>5; pos2(4) = pos(4)+5; end   
        A(pos2(2):pos2(4),pos2(1):pos2(3))=0;
    end
    end
    CountImage=CountImage+A;
    Combined2=Combined2+imdilate(A,strel('disk',2));
end


Combined2(Combined2<1)=+Inf;  %num minimal count to keep roi
X=((Combined)./Combined2);
X(X<=1.5)=0;
figure;
X=-X;
imshow(-X,[0 4])
X(X==0)=-Inf;

Y=X;


X = watershed(X)-1;
X=imdilate(X,strel('disk',1));


for(i=2:size(X,1)-1) 
    for(j=2:size(X,2)-1) 
        ref=X(i,j);
        if(ref~=0)
            if(((X(i,j+1)~=ref)&&(X(i,j+1)~=0))||((X(i+1,j)~=ref)&&(X(i+1,j)~=0))||((X(i,j-1)~=ref)&&(X(i,j-1)~=0))||((X(i-1,j)~=ref)&&(X(i-1,j)~=0)))
                X(i,j)=0;
            end
        end
    end
end

totcell=max(max(X));
removed=0;
CellAreas=[];
Masks=zeros(size(Images(:,:,1)));
for(counter=1:totcell)
   if((size(find(X==counter),1)>5 && size(find(X==counter),1)<200)&&~min(Y(X==counter)==-inf))
       CellAreas(counter-removed)=size(find(X==counter),1);
       Masks(X==counter)=counter-removed;
   else
       removed=removed+1;
   end
end

rgb = label2rgb(Masks,'jet',[.5 .5 .5]);
MasksDilated=imdilate(Masks,strel('disk',5));
Total=0;

for(i=1:max(max(Masks)))
    Inter=imdilate((Masks==i),strel('disk',14))-imdilate((Masks==i),strel('disk',6));
    Inter(MasksDilated~=0)=0;
    NeuropilMasks(:,:,i)=Inter;
    CellMasks(:,:,i)=double((Masks==i));
end

for(i=1:max(max(Masks)))
    Total=Total+NeuropilMasks(:,:,i);
end

figure;
subplot(2,2,1);
imshow(mean(Images,3),[0 prctile(prctile(mean(Images,3),90),90)]);
subplot(2,2,2);
imshow(CountImage,[0 1]);
subplot(2,2,3);
imshow(rgb);
subplot(2,2,4);
imshow(Total+Masks,[0 1]);

%%
icaguidata = struct('ica', [], 'movm', meanim, 'movcorr', ln_meanim, 'snrsort', []);

icaguidata.pars.threshold = Threshold;
icaguidata.pars.n = n; icaguidata.pars.m = m;

icaguidata.pars.downsample_xy = p.downsample_xy;

% Pull ROI signal
signal = sbxALSimplePullSignalsFromMovie(mov,logical(CellMasks),p.downsample_xy); %% Pulls signals from a movie, looping over nrois, without chunking movie
neuropilsignal = sbxALSimplePullSignalsFromMovie(mov,logical(NeuropilMasks),p.downsample_xy); %% Pulls signals from a movie, looping over nrois, without chunking movie


for i = 1:size(CellMasks,3)
    icaguidata.ica(i).filter = meanim ./ CellMasks(:,:,i);
    icaguidata.ica(i).filter(isinf(icaguidata.ica(i).filter)) = 0;
    icaguidata.ica(i).filter(isnan(icaguidata.ica(i).filter)) = 0;
    icaguidata.ica(i).trace = signal(:,i);
    icaguidata.ica(i).neuropiltrace = neuropilsignal(:,i);
end

save(icapath, 'icaguidata', '-v7.3');
end