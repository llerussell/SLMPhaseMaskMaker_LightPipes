function SLMPhaseMaskMaker(varargin)
% Adam Packer December 2 2013, Lloyd Russell April 2014, Jan 2015
% varargin can be:
%   none, an open file dialog box will open
%   path to image file on disk
%   loaded image data, save name must be provided as second argument

% Load targets image
if nargin == 0  % open file dialog if no filepath provided
    [FileName,PathName] = uigetfile('*.tif*','Select the targets file');
    filepath = [PathName filesep FileName];
    cd(PathName)
    targetsImg = imread(filepath);
    transform_choice = questdlg('Transform points?', 'User input required', 'Yes','No', 'Yes');
elseif ischar(varargin{1})  % filepath provided
    filepath = varargin{1};
    targetsImg = imread(filepath);
    transform_choice = 'Yes';
else  % assume data array provided first, and savename second
    targetsImg = varargin{1};
    filepath = varargin{2};
    transform_choice = 'Yes';
end
[~,filename] = fileparts(filepath);

% Find target coordinates
[y,x,~] = find(targetsImg);

% Load previously calibrated transform
load('20150421_tform_2P_SLM')

% Convert target *coordinates* directly from 2P to SLM space
switch transform_choice
    case 'Yes'
        [u,v] = transformPointsForward(tform,x,y);
        u     = round(u);
        v     = round(v);
    case 'No'
        u = x; v = y;
end

% Build transformed targets image
SLMtargets = uint8(zeros(512,512));
for i = 1:length(u)
    SLMtargets(v(i),u(i)) = 255;
end
imwrite(SLMtargets, [filename '_Transformed.tif']);

% Weight the target pixel intensity by distance from zero order
slope                = 75/183;
b                    = 25;
DistFromZeroOrderIdx = 1:0.1:183;
WeightByDist         = slope*DistFromZeroOrderIdx + b;
WeightByDist         = WeightByDist/100;
for i = 1:length(u)
    d(i) = sqrt((u(i)-256)^2 + (v(i)-256)^2);
end
for i = 1:length(u)
    ThisDistIdx = find(d(i)<DistFromZeroOrderIdx,1,'first');
    ThisDistIdx(d(i)>max(DistFromZeroOrderIdx)) = length(DistFromZeroOrderIdx); %added 20141129
    p(i) = 255*WeightByDist(ThisDistIdx);
end

% Build final transformed, weighted targets image
SLMtargets = uint8(zeros(512,512));
for i = 1:length(u)
    SLMtargets(v(i),u(i)) = p(i);
end
imwrite(SLMtargets, [filename '_TransformedWeighted.tif']);

% Computer generated hologram using the Gerchberg-Saxton algorithm
% Dr F.A. van Goor, University of Twente. April 2010
% Reorganized by Adam Packer December 2013
Nitt    = 10; % number of itererations
m       = 1;
nm      = 1e-9*m;
mm      = 1e-3*m;
cm      = 1e-2*m;
lambda  = 1064*nm; %red HeNe laser
SLMsize = 7.68*mm; %the HoloEye LCD is W x H = 25.4mm x 19.05 mm
N       = 512; %The HoloEye LCD has W x H = 800 x 600 pixels, we use a square grid
imageD  = double(SLMtargets); %The LightPipes command 'LPSubIntensity' requires an array of doubles as input
UniformIntensity = ones(N); % We need a matrix filled with 1's to substitute a uniform intensity profile

% The iteration loop to get the phase distribution
F = LPBegin(SLMsize,lambda,N); % We start with a plane wave distribution with amplitude 1 and phase 0
for i = 1:Nitt
    F = LPPipFFT(1,F); %Take the 2-D Fourier transform of the field
    F = LPSubIntensity(imageD,F); % Substitute the original intensity distribution, leave the phase unchanged
    F = LPPipFFT(-1,F); % Take the inverse Fourier transform
    F = LPSubIntensity(UniformIntensity,F); % Substitute a uniform intensity, leave the phase unchanged
    fprintf('%s ','.'); % monitor the number of iterations done in the command window
end;
fprintf('\n');

Phase       = LPPhase(F); %Extract the phase distribution from the field
PhaseZeroed = Phase+abs(min(min(Phase)));

% Convert
% Phase8      = PhaseZeroed*(255/max(max(PhaseZeroed)));
% phaseMask8  = uint8(Phase8);
Phase16     = PhaseZeroed*(65535/max(max(PhaseZeroed)));
phaseMask16 = uint16(Phase16);

% Save
imwrite(phaseMask16, [filename '_LPphase.tiff']);
