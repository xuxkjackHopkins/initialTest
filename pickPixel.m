function pixel_index = pickPixel(image,threshold)
% This function is used to pick pixel index of bioluminescence signal in 2D
% BLIs.
% input : 
%    image : 2D BLI image
%    background: background of the image.
%    threshold : pixel value threshold to binarize BLI image
% output:
%    pixel_index : index of pixel (u, v)
% Jack Xu, Jan 2020, Johns Hopkins University.
bw = image >= threshold; % treshold used to select pixels
bw_filt = bwareafilt(bw,[4, size(bw,1)*size(bw,2)]); % remove pixels less than 4
[u, v] = find(bw_filt>0);
pixel_index = [u v];
