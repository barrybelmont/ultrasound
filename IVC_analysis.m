
%% Open dataset
dataset = 'dialysis';
cd(strcat('/Volumes/usb/',dataset))
files = dir;
files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'),files));

for i = 1:numel(files)
    % Patient ID
    data(i).name = files(i).name(1:7);
    
    % Look for hyphenated information (before/after)
    if length(find(files(i).name=='-',2)) > 1
        hyphens = find(files(i).name=='-',2);
        data(i).when = files(i).name(hyphens(1)+1:hyphens(2)-1);
    end
        
end

%% Go through samples
% Which file number are we working on now?
sample = 10;

% Read in DICOM, get size/sampling information
dicomFile = dicomread(files(sample).name);
[height, width, ~, nFrames] = size(dicomFile);
info = dicominfo(files(sample).name);
MM_PER_PIXEL = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX*10;
dtUS = (info.FrameTime)*0.001;
FsUS = 1/dtUS;
time = (1:nFrames)/FsUS;

% Display the (first) image and pick points to track
indFrame = 1;
objectFrame = dicomFile(:,:,:,indFrame);
imshow(objectFrame)
title('Select points along the vessel edge, then hit "Enter"')
figHandle = gcf;
[poiX, poiY] = getpts(figHandle);
close

poiX = round(poiX);     poiY = round(poiY);
points = [poiX, poiY];
nPoints = size(poiX,1);
pointLog = zeros(nPoints, 2, nFrames);

% Track points of interest
tracker = vision.PointTracker('MaxBidirectionalError', inf);
initialize(tracker, points(:,:,1), objectFrame);

while indFrame <= nFrames
    % Track the points
    frame = dicomFile(:,:,:,indFrame);
    [points, validity] = step(tracker, frame);
    pointLog(:,:,indFrame) = points;
    
    indFrame = indFrame + 1;
    
end

% Diameter over time
pointDist = permute(sqrt((pointLog(1:2:end,1,:) - pointLog(2:2:end,1,:)).^2 ...
    + (pointLog(1:2:end,2,:) - pointLog(2:2:end,2,:)).^2), [3 1 2]);


% Color for points
color = [0,    0.4470,    0.7410;
    0.8500,    0.3250,    0.0980;
    0.9290,    0.6940,    0.1250;
    0.4940,    0.1840,    0.5560;
    0.4660,    0.6740,    0.1880;
    0.3010,    0.7450,    0.9330;
    0.6350,    0.0780,    0.1840];

% Include a reference image with points tracked
out = insertMarker(repmat((rgb2gray(objectFrame(:,:,:,1))), [1 1 3]), ...
    [poiX(1:2), poiY(1:2)], 'o', 'Color', color(1,:)*255, 'size', 15);
for np = 3:2:nPoints
    out = insertMarker(out, [poiX(np:np+1, :), poiY(np:np+1, :)],...
        'o', 'Color', color((np-1)/2+1,:)*255,'size',15);
end

data(sample).image = out;

% Create data structure to report the point pairs, the diameter over time,
% the maximum diameter, the minimum diameter, and the overall dIVC
for ind = 1:floor(nPoints/2)
    data(sample).(strcat('pair',num2str(ind))) = ...
        [poiX(2*ind-1), poiY(2*ind-1), poiX(2*ind), poiY(2*ind)];
    data(sample).(strcat('diam',num2str(ind))) = ...
        pointDist(:,ind)*MM_PER_PIXEL;
    data(sample).(strcat('max',num2str(ind))) = ...
        max(pointDist(:,ind)*MM_PER_PIXEL);
    data(sample).(strcat('min',num2str(ind))) = ...
        min(pointDist(:,ind)*MM_PER_PIXEL);
    data(sample).(strcat('dIVC',num2str(ind))) = ...
        (max(pointDist(:,ind)*MM_PER_PIXEL) ...
        - min(pointDist(:,ind)*MM_PER_PIXEL))...
        ./max(pointDist(:,ind)*MM_PER_PIXEL)*100;   % In percent
end

data(sample).diamALL = pointDist*MM_PER_PIXEL;
data(sample).MM_PER_PIXEL = MM_PER_PIXEL;
