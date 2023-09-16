% INSTRUCTIONS:
%       Modify datapath, fname, loadROI, slice
% OUTPUT:
%       T1_mean, T1_pixel and masks of T1
clear all;
clc; 
close all;

addpath('toolbox');

datapath = ['..' filesep 'exampledata'];

fname = [datapath,filesep,'DFA.mat']
load(fname)


TR=25/1e3;  %   TR = 25 ms


a5 = 5/180*pi;
S5_a_r = img5./a5;
S5_a_p = img5.*a5;



a30 = 30/180*pi;
S30_a_r = img30./a30;
S30_a_p = img30.*a30;
T1map = 2.*TR.*(S5_a_r - S30_a_r)./(S30_a_p - S5_a_p);  %   T1 in the unit of s

%%
loadROI = 1; % change to 1 when you have plotted the ROI!!
ROI1 = [datapath filesep 'MASK.mat'];

% DRAW or INPUT all the ROI (totally 15 slices here)
nslice = size(T1map, 3);
if loadROI == 1
    load(ROI1);
else
    for i = 1:nslice
        displayimg=squeeze(T1map(:,:,i));
        scrsz = get(0,'ScreenSize');
        h2=figure;
        set(h2,'Position',[1 scrsz(4)*0.6 scrsz(3)*0.4 scrsz(4)*0.3]);
        imagesc(displayimg, [0 2]);
        colorbar('location','Eastoutside','FontSize', 18)
        colormap(gray(255))
        set(gca,'dataAspectRatio', [1 1 1]);
        axis image;
        axis off;
        ROI(:,:, i)=logical(roipoly);
        tmpR=displayimg.*ROI(:,:,i);
        ROI(:,:,i)=((1<tmpR)&(tmpR<1.6));
        
    end
        mask = ROI;
        save(ROI1,'mask');
end

if exist('mask','var') == 1
    ROI = mask;
end

for i = 1:nslice
    displayimg = squeeze(T1map(:,:,i));
    single_mask = ROI(:, :, i);
   
    % T1_mean for each slice
    T1mean(i) = mean(displayimg(logical(single_mask)));
end

T1mappixel = zeros(size(ROI,1), size(ROI,2), nslice);
for i = 1:nslice
    for idx = 1:size(ROI,1)
        for idy = 1:size(ROI, 2)
            if ROI(idx, idy, i) == 1
       T1mappixel(idx, idy, i) = T1map(idx,idy,i);
            end
        end
    end
end

sname = 'T1mappixel.mat';
sname = [datapath filesep sname];
save(sname, 'T1mappixel');
