function [newPoints] = findNearestHarrisCorner(I, points, thresh)
%FINDNEARESTHARRISCORNER Finds the nearest Harris corner satisfying a
%specified threshold to a selected point
%   Input:
%   I           - the image in which to find corners
%   points      - the points about which to look for points
%   thresh      - the threshold to determine a Harris corner
%
%   Output:
%   newPoints   - the new points to track


% Determine if it is a grayscale or RGB image and convert to grayscale
% using doubles
if size(I,3) == 3
    I = double(rgb2gray(I));
else I = double(I);
end

% Perform convolutions for Harris corners
blur    = [1 6 15 20 15 6 1];
blur    = blur / sum(blur);
prefilt = [0.223755 0.552490 0.223755];
derivfilt = [-0.453014 0 0.45301];
fx     = conv2( conv2( I, prefilt', 'same' ), derivfilt, 'same' );
fy     = conv2( conv2( I, prefilt, 'same' ), derivfilt', 'same' );
fx2    = conv2( conv2( fx .* fx, blur', 'same' ), blur, 'same' );
fy2    = conv2( conv2( fy .* fy, blur', 'same' ), blur, 'same' );
fxy    = conv2( conv2( fx .* fy, blur', 'same' ), blur, 'same' );

m = (fx2 + fy2)/2;
d = fx2.*fy2 - fxy.^ 2;
n = sqrt(m.^2 - d);

% Find the 'quality' of the Harris corner
quality = min(abs(m - n),abs(m + n));

% Search for closest best Harris corner
for indPoints = 1:size(points,1)
    maxVal = 0;
    a = 1;
    while maxVal < (thresh - a) % Drop the quality as you search further (this needs to be improved)
        A = quality((points(indPoints,2)-a):(points(indPoints,2)+a),(points(indPoints,1)-a):(points(indPoints,1)+a));
        
        [maxVal, maxInd] = max(A(:));
        [x, y] = ind2sub(size(A),maxInd);
        
        origPoint = ((((2*a)+1)^2-1)/2)+1;
        [xc, yc] = ind2sub(size(A),origPoint);
        
        a = a+1;
    end
    newPoints(indPoints,:) = [points(indPoints,1)+(x-xc), points(indPoints,2)+y-yc];
    
end

end



