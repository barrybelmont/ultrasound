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
            yoi = round(yfit+rfit+5/MM_PER_PIXEL);
            
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
filename = 'C0029600';
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

clear ypts xpts
for indCol = round(2*width/8):round(6*width/8)
    if size(find(Iseg1(round(height/2):height, indCol) == 2, 1, 'first'),1) == 1
        ypts(indCol) = round(height/2) + find(Iseg1(round(height/2):height, indCol) == 2, 1, 'first');
        xpts(indCol) = indCol;
    end
end

[xfit, yfit, rfit] = circfit(xpts,ypts);
imagesc(Iseg1)
hold on
rectangle('position',[xfit-rfit,yfit-rfit,rfit*2,rfit*2],...
    'curvature',[1,1],'linestyle','-','edgecolor','y','LineWidth',2);


%% Fit circle
[xfit, yfit, rfit] = circfit(xpts,ypts);

xoi = round(xfit);
yoi = round(yfit+rfit+3/MM_PER_PIXEL);

if Iseg1(yoi,xoi) == 2
    for ind = yoi:height-10
        if size(find(Iseg2(ind, 1:xoi) == 1, 2, 'last'),1) == 1
            edges(ind) = find(Iseg2(ind, 1:xoi) == 1, 1, 'last');
        end
    end
end

leftWall = find(Iseg2(yoi, 1:xoi) == 3, 1, 'last');
rightWall = round(xoi+find(Iseg2(yoi, xoi:end) == 3, 1, 'first'));

ONSD = (rightWall - leftWall)*MM_PER_PIXEL;

%%

RGB = insertMarker(imageROI,[leftWall yoi; rightWall yoi],'o','Color','yellow','size',8);
imshow(RGB)
hold on

plot([leftWall, rightWall], [yoi,yoi],'y-','LineWidth',1)


%%
clear Iseg RGB          % changed dimensions from those above
thresh = multithresh(I,threshLevel(1));
Iseg(:,:) = imquantize(I,thresh);
RGB(:,:,:) = label2rgb(Iseg);

[height, width] = size(Iseg);

for indCol = 20:200
    ypts(indCol) = round(height/2) + find(Iseg(round(height/2):height, indCol) == 2, 1, 'first');
end


[xfit, yfit, rfit] = circfit(20:200, ypts(20:200));

figure
plot(20:200,ypts(20:200))
hold on
rectangle('position',[xfit-rfit,yfit-rfit,rfit*2,rfit*2],...
    'curvature',[1,1],'linestyle','-','edgecolor','r');
% plot(xfit,yfit,'r')


%%



info = dicominfo(filename);
MM_PER_PIXEL = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY*10;

xoi = xfit;
yoi = round(yfit+rfit+5/MM_PER_PIXEL);

thresh = multithresh(I,threshLevel(2));
Iseg = imquantize(I,thresh);

find(Iseg(yoi, 1:xoi) == 3, 1, 'last')
round(xoi+find(Iseg(yoi, xoi:end) == 3, 1, 'first'))



%%

imshow(imageROI)
hold on
rectangle('position',[xfit-rfit,yfit-rfit,rfit*2,rfit*2],...
    'curvature',[1,1],'linestyle','-','edgecolor','y','LineWidth',2);

