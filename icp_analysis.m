%% RUN THROUGH ALL DATA
dataset = 'icp';
cd('/Volumes/usb/')


cd(strcat('/Volumes/usb/',dataset))
folders = dir;
folders = folders(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders));


counter = 1;
tic
for i = 1:numel(folders)
    cd(strcat('/Volumes/usb/',dataset,'/',folders(i).name))
    files = dir;
    files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'),files));
    
    for j = 1:numel(files)
        filename = files(j).name;
        info = dicominfo(filename);
        if isfield(info,'SequenceOfUltrasoundRegions')
            MM_PER_PIXEL = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY*10;
            
            icp(counter).PATIENT_ID = folders(i).name;
            icp(counter).filename = files(j).name;
            
            images = dicomread(filename);
            
            ROI = [info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0, ...
                info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1, ...
                info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0, ...
                info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1];
            
            imageROI = rgb2gray(images(ROI(3):ROI(4),ROI(1):ROI(2),:,1));
            
            Iseg1 = imquantize(imageROI,multithresh(imageROI,1));
            Iseg2 = imquantize(imageROI,multithresh(imageROI,2));
            
            [height, width] = size(Iseg1);
            
            for indCol = round(width/4):round(3*width/4)
                if size(find(Iseg1(round(height/2):height, indCol) == 2, 1, 'first'),1) == 1
                    ypts(indCol) = round(height/2) + find(Iseg1(round(height/2):height, indCol) == 2, 1, 'first');
                    xpts(indCol) = indCol;
                end
            end
            
            [xfit, yfit, rfit] = circfit(xpts,ypts);
            
            xoi = round(xfit);
            yoi = round(yfit+rfit+3/MM_PER_PIXEL);
            
            leftWall = find(Iseg2(yoi, 1:xoi) == 3, 1, 'last');
            rightWall = round(xoi+find(Iseg2(yoi, xoi:end) == 3, 1, 'first'));
            
            ONSD = (rightWall - leftWall)*MM_PER_PIXEL;
            
            icp(counter).ONSD = ONSD;
            icp(counter).xoi = xoi;
            icp(counter).yoi = yoi;
            icp(counter).leftWall = leftWall;
            icp(counter).rightWall = rightWall;
            icp(counter).MM_PER_PIXEL = MM_PER_PIXEL;
            
            counter = counter + 1;
            %
            
        end
    end
    
end
toc

%% DEBUGGING CODE STARTS FROM HERE ON DOWN
filename = 'C0028341';
images = dicomread(filename);
info = dicominfo(filename);
MM_PER_PIXEL = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY*10;

ROI = [info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0, ...
    info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1, ...
    info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0, ...
    info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1];

imageROI = rgb2gray(images(ROI(3):ROI(4),ROI(1):ROI(2),:,1));
imshow(imageROI)

%% Create segmented images
Iseg1 = imquantize(imageROI,multithresh(imageROI,1));
Iseg2 = imquantize(imageROI,multithresh(imageROI,2));

%% Look for bottom edge
[height, width] = size(Iseg1);

clear pts
buffer = 20;

for indCol = 1:width
    if size(find(imageROI(round(height/2):height-buffer,indCol) >= 20, 1, 'first'),1)==1
        pts(indCol) = round(height/2) + ...
            find(imageROI(round(height/2):height-buffer,indCol) >= 20,1, 'first');
    end
end

[ylo, xlo] = max(pts);

yoi = ylo + 3/MM_PER_PIXEL;
xoi = xlo;


