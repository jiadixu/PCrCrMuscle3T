% INSTRUCTIONS:
%       1. Modify datapath, fname, ROI1, xname;
%       2. Modify B1 value (pwr), saturation time (tsat)
%       3. LOAD pixel-wise T1 map 
%
% OUTPUT: 
%       1. Pixel-wise CEST Z-spectrum
%       2. Amide, PCr and Cr map for each slice 

clear all;
clc; 
close all;

addpath(['..' filesep 'code' filesep 'toolbox']);

%% what you may need to modify
datapath = ['..' filesep 'exampledata'];

tube = 1; % how many tubes do you want to extract at one time

fname = 'CEST2s.mat';
fname=[datapath filesep fname];
%load CEST images
load(fname);
%load offset list
offsetname=[datapath filesep 'crlist.mat']; % ensure the parameter is named as FreqHz
load(offsetname);
% load MASK 
loadROI = 1; % change to 1 when you have plotted the ROI!
ROI1 = [ datapath filesep 'MASK.mat']; 

%%load T1map
tname = [datapath filesep 'T1mappixel'];
load(tname);

%%save directory
savedatadir=[datapath filesep 'FitOutput']
if exist(savedatadir,'file')==0
    mkdir(savedatadir);
end


IMGtmp = squeeze(imgs);
IMGtmp=imresize(IMGtmp, [64 64]);


if loadROI == 1
    load(ROI1);
else
    for i = 1: tube
        ROI(:,:,i)=roipoly;
        tmpR=displayimg.*ROI(:,:,i);
        ROI(:,:,i)=(tmpR>6e+06);
    end
    save(ROI1,'ROI');
end


amidemap = zeros(size(IMGtmp,1),size(IMGtmp,2),size(IMGtmp,3)-2);
Crmap = zeros(size(IMGtmp,1),size(IMGtmp,2),size(IMGtmp,3)-2);
PCrmap = zeros(size(IMGtmp,1),size(IMGtmp,2),size(IMGtmp,3)-2);
for itube =3:size(IMGtmp, 3)-1 %only plot the ten slices in the middle
    for idx = 1:size(IMGtmp,1)
        for idy = 1:size(IMGtmp,2)
               if mask(idx, idy, itube+1)
                    CEST(:, 1) = squeeze(IMGtmp(idx,idy,itube,:));
                    CESTM0(:, 2)=CEST(:,1)./CEST(2,1);
                    CESTM0(:, 1) = FreqHz/128;
                    [M, d] = min(CESTM0(:, 2));
                    mini = CESTM0(d,1);
                    minamp = CESTM0(d,2);
                    for k=1:size(CESTM0,1)
                        CESTM0(k,1)=CESTM0(k,1)-mini;
                        CESTM0(k,2)=(CESTM0(k,2)-minamp)/(1-minamp);
                    end
                    CESTpixel(idx, idy, itube-2,:, :) = CESTM0;
                    Freq = CESTM0(:,1);
                    Z_spectrum = CESTM0(:,2);
                    T1mean = T1mappixel(idx, idy, itube+1);
                    FitParam.R1 = 1/T1mean;
                    FitParam.WholeRange = [-1.5, 1.5];
                    FitParam.ifshowimage = 0; % If =1, show Z-spectrum
                    B0 = WASSR(Freq, Z_spectrum, FitParam);
                    Freq = CESTM0(:,1)-B0;

                    %interpolate 10X data points in CEST
                    xraw=CESTM0(3:end, 1)';
                    yraw=CESTM0(3:end, 2)';
                    newx=-2:0.01:8;
                    newy=interp1(xraw,yraw,newx);
                    for i=1:size(newx,2)/10
                    CESTnew(i,1)=newx(1,i*10);
                    CESTnew(i,2)=newy(1,i*10);
                    end
                    Freq=CESTnew(:,1);
                    Z_spectrum=CESTnew(:,2);

            
                    FitParam.WholeRange = [0.8, 7.8];    % total fitting range
                    FitParam.PeakOffset = 2; %this define the background offset, any value between 1.5-5 ppm will be fine 
                    FitParam.PeakRange = [1.5, 5]; % fitting CEST peaks range
                    FitParam.Magfield = 42.58*3; % Lamor frequency for 3T
                    FitParam.Offset1=[3.5 3.5]; % peak 1 offset range
                    FitParam.Offset2=[2 2]; % peak 2 offset range
                    FitParam.Offset3=[2.5 2.7]; % peak 3 offset range

                    FitParam.satpwr = 0.6; % saturation power (uT)
                    FitParam.tsat = 2; %saturation time (s)

                    [FitResult, FitParam] = PLOFpixel(Freq, Z_spectrum, FitParam);

                    amidemap(idx, idy, itube-2) = FitResult.DeltaZpeak1;
                    Crmap(idx, idy, itube-2) = FitResult.DeltaZpeak2;
                    PCrmap(idx, idy, itube-2) = FitResult.DeltaZpeak3;
                                   
               end
               
          end
    end
end
%% 
% Section 0: Plot PCr abs map
PCrabsmap=PCrmap.*(39.58/0.026); %calibrated by 31P MRS results
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
montage(PCrabsmap, 'DisplayRange',[0,45],'Size',[2, 5]);
colorbar('location','Eastoutside','FontSize', 18);
title('PCr Concentration map (mM)');
set(gca,'dataAspectRatio', [1 1 1]);
colormap(inferno);
axis image;
axis off;

xname = [savedatadir filesep  'PCrConMap']; % save name 
saveas(gcf,[xname,'.tif'],'tiff'); 

% Section 1: Plot the Crmap
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
montage(Crmap, 'DisplayRange',[0,0.03],'Size',[2, 5]);
colorbar('location','Eastoutside','FontSize', 18);
title('CrCEST map');
set(gca,'dataAspectRatio', [1 1 1]);
colormap(inferno);
axis image;
axis off;

xname = [savedatadir filesep  'CrCESTMap']; % save name 
saveas(gcf,[xname,'.tif'],'tiff'); 

%Plot the PCrmap
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
montage(PCrmap, 'DisplayRange',[0,0.03],'Size',[2, 5]);
colorbar('location','Eastoutside','FontSize', 18);
title('PCrCEST map');
set(gca,'dataAspectRatio', [1 1 1]);
colormap(inferno);
axis image;
axis off;

xname = [savedatadir filesep  'PCrCESTMap']; % save name 
saveas(gcf,[xname,'.tif'],'tiff'); 

%Plot the amidemap
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
montage(amidemap, 'DisplayRange',[0,0.03],'Size',[2, 5]);
colorbar('location','Eastoutside','FontSize', 18);
title('AmideCEST map');
set(gca,'dataAspectRatio', [1 1 1]);
colormap(inferno);
axis image;
axis off;

xname = [savedatadir filesep  'AmideCESTMap']; % save name 
saveas(gcf,[xname,'.tif'],'tiff'); 

dname = [ savedatadir filesep '3DEPI2s06uT.mat']; % save name
save(dname, 'CESTpixel', 'FitParam');