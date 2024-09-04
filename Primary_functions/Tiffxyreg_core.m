function shifts = Tiffxyreg_core(im, ref, gausssize)
%Core function for Tiffxyreg
% shifts = Tiffxyreg_core(im, ref, gausssize)

% Initialize
im2 = zeros(size(im));

% Local normalize
for i = 1 : size(im, 3)
    frame = medfilt2(im(:,:,i), [2 2], 'symmetric');
    f_prime = frame - imgaussfilt(frame, gausssize(1));
    frame = f_prime ./ (imgaussfilt(f_prime.^2, gausssize(2)).^(1/2));
    frame(isnan(frame)) = 0;
    im2(:,:,i) = frame;    
end

% Get shifts
[shifts,~] = stackRegisterMA_RR(im2,ref);
    
end