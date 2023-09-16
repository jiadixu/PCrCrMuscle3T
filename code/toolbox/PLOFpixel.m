function [FitResult, FitParam] = PLOFpixel(Offset, Saturation, FitParam)
% Fit the Z-spectrum using polynomial and Lorentzian line-shape fitting (PLOF) method
%
% If you have used this function in a scientific publication,we would
% appreciate citations to the following papers:

% [1]	Chen L, Wei Z, Cai S, Li Y, Liu G, Lu H, Weiss RG, van Zijl PCM, Xu
% J. High-resolution creatine mapping of mouse brain at 11.7 T using
% non-steady-state chemical exchange saturation transfer. NMR Biomed
% 2019:e4168. 
% [2]	Chen L, Barker PB, Weiss RG, van Zijl PCM, Xu J.
% Creatine and phosphocreatine mapping of mouse skeletal muscle by a
% polynomial and Lorentzian line-shape fitting CEST method. Magn Reson Med
% 2019;81(1):69-78. 
% [3]	Chen L, Zeng H, Xu X, Yadav NN, Cai S, Puts NA,
% Barker PB, Li T, Weiss RG, van Zijl PCM, Xu J. Investigation of the
% contribution of total creatine to the CEST Z-spectrum of brain using a
% knockout mouse model. NMR Biomed 2017;30(12):e3834.
%
% Please contact Lin Chen at chenlin@jhu.edu or  chenlin0430@gmail.com if you have any questions about the code. 

% 1/10/2022 UPDATED:
%       add FitResult.SNR
% 1/21/2022 UPDATED:
%       modify FitResult.SNR as mean-of-the-middle / std-of-the-squares-on-four-corners
% 2/4/2022 UPDATED:
%       modify FitResult.SNR as to shrink the noise mask from 19*19 to
%       16*16, since 16*16*4 = 1024, still fulfills the requirement of 1000
% 2/15/2022 UPDATED:
%       add goodnessOfFit for background fitting
% 2/21/2022 UPDATED:
%       change Offset(idx) to FitResult.xindex(idx)
% 8/21/2023 UPDATED:
%       Jiadi Xu add three peaks and control three peak fitting range

warning off
if size(Offset,1) == 1
    Offset = Offset';
end

if size(Saturation, 1) == 1
    Saturation = Saturation';
end

if (Saturation(1)) > 10
    error('Z-spectrum have to be 0-1, can not be percentage')
end


[Offset, Saturation] = Extract_Zspectrum(Offset, Saturation, FitParam.WholeRange);
[Offset_background, Saturation_background] = CutOff_Zspectrum(Offset, Saturation, FitParam.PeakRange);

FitResult.Offset = Offset;
FitResult.Saturation = Saturation;
FitResult.Offset_background = Offset_background;
FitResult.Saturation_background = Saturation_background;

%-----------------------------------------------%
%-------------- Two-step fitting ---------------%
%-----------------------------------------------%

FitResult.xindex = linspace(min(Offset), max(Offset), 100)';

x0 = [0.1, 1, 0.5, -19, 1,... %background functions
        0.1, 1, FitParam.Offset1(1), ...     % amide peak fix at 3.5 ppm
        0.1, 1, FitParam.Offset2(1), ...       % Guan peak fix at 2 ppm
        0.1, 1, FitParam.Offset3(1)];       % PCr or amine peak fix at 2.5-2.7 ppm
lb = [0, 0, 0, -100, -10,...
        1e-8, 0.4, FitParam.Offset1(1), ...
        1e-8, 0.4, FitParam.Offset2(1), ...
        1e-8, 0.6, FitParam.Offset3(1)]; 
ub = [200, 100, 100, 100, 10,...
            10, 1, FitParam.Offset1(2), ...
            10, 1, FitParam.Offset2(2), ...
            10, 2, FitParam.Offset3(2)]; 
        
%-------------background fitting --------------%

FitParam.peak = 0;

options = optimset('MaxFunEvals', 1e6,'TolFun', 1e-6,'TolX', 1e-6, 'Display',  'off' );
% f = @(x, Offset_background)CurveFunction(x, Offset_background, FitParam);
% [FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset_background, Saturation_background, lb, ub, options, FitParam);
% [FitResult.Coefficents, resnorm] = lsqcurvefit(f, x0, Offset_background, Saturation_background, lb, ub, options);
[FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset_background, Saturation_background, lb, ub, options, FitParam);

% fitting result of Z-spectrum and R-spectrum
[FitResult.Background, Rbak] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);
Backgroundtmp = CurveFunction(FitResult.Coefficents, Offset, FitParam);
FitResult.DeltaZspectrum = Backgroundtmp - Saturation; 


% R^2 for the evaluation of fitting
FitResult.RsquareBG = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset_background, FitParam), Saturation_background, 'NMSE');

%-------------CEST all-peak fitting new--------------%
FitParam.peak = 4;

% fix the background parameters
lb(1: 5) = FitResult.Coefficents(1 : 5); 
ub(1: 5) = FitResult.Coefficents(1 : 5); 

[FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
[FitResult.Curve, Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);

FitResult.DeletaFitZ = FitResult.Background-FitResult.Curve;
FitResult.DeletaFitR = Rtmp - Rbak;

% %-------------CEST all-peak fitting old --------------%
% FitParam.peak = 1;
% 
% % fix the background parameters
% lb(1: 5) = FitResult.Coefficents(1 : 5); 
% ub(1: 5) = FitResult.Coefficents(1 : 5); 
% 
% [FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
% [FitResult.Curve, Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);
% 
% FitResult.DeletaFitZ = FitResult.Background-FitResult.Curve;
% FitResult.DeletaFitR = Rtmp - Rbak;

%------------ Calculate deltaZ and assign the Coefficents ---------------%

FitResult.Rpeak1 = FitResult.Coefficents(6);% R_exch
FitResult.FitPeakOffset1 = FitResult.Coefficents(8); % delta_w_exch
FitResult.Rpeak2 = FitResult.Coefficents(9);% R_exch
FitResult.FitPeakOffset2 = FitResult.Coefficents(11); % delta_w_exch
FitResult.Rpeak3 = FitResult.Coefficents(12);
FitResult.FitPeakOffset3 = FitResult.Coefficents(14);

FitResult.MT = FitResult.Coefficents(1); % C_0
FitResult.RsquareAll = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset, FitParam), Saturation, 'NMSE');


%------------ Calculate DeltaZ-spectrum ---------------%

% amide peak
FitParam.peak = 2;

[Ztmp Rtmp]= CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZamideFit =FitResult.Background -Ztmp;
[FitResult.DeltaZpeak1,  idx]= max(FitResult.ZamideFit);
FitResult.DeltaZpeak1Offset = FitResult.xindex(idx); % Offset(idx);
FitResult.RamideFit =Rtmp - Rbak;
[FitResult.DeltaRpeak1, idx] = max(FitResult.RamideFit);
FitResult.DeltaRpeak1Offset = FitResult.xindex(idx); % Offset(idx);

% guan peak
FitParam.peak = 3;

[Ztmp Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZguanFit =FitResult.Background -Ztmp;
[FitResult.DeltaZpeak2,  idx] = max(FitResult.ZguanFit);
FitResult.DeltaZpeak2Offset = FitResult.xindex(idx); % Offset(idx);
FitResult.RguanFit =Rtmp - Rbak;
[FitResult.DeltaRpeak2,  idx] = max(FitResult.RguanFit);
FitResult.DeltaRpeak2Offset = FitResult.xindex(idx); % Offset(idx);

% amine peak
FitParam.peak = 5;

[Ztmp Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZamineFit =FitResult.Background -Ztmp;
[FitResult.DeltaZpeak3,  idx] = max(FitResult.ZamineFit);
FitResult.DeltaZpeak3Offset = FitResult.xindex(idx); % Offset(idx);
FitResult.RamineFit =Rtmp - Rbak;
[FitResult.DeltaRpeak3,  idx] = max(FitResult.RamineFit);
FitResult.DeltaRpeak3Offset = FitResult.xindex(idx); % Offset(idx);




   if FitParam.ifshowimage == 1
     PlotFitResult(FitResult,FitParam);
        end
end