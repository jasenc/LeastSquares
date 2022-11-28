clear all
close all
clc

% reads in the laser range finder data 
% the first column of the data denotes the angle
% the second column of the data denotes the range or distance in mm
data = load('lineData.txt');

% converts the data in polar coordinates to X and Y coordinates
Xdata = data(:,2).*cos(data(:,1));
Ydata = data(:,2).*sin(data(:,1));

% plots the data
plot(Xdata, Ydata, 'r+')
axis([-1000 1500 -1000 1500])
hold on
%% We will first separate the data set into M different clusters.  We will
%% then try to fit lines to each individual cluster

tol = 100;  % if two consecutive points are more than 100 mm (10 cm) apart, start a new cluster
K = 1;      % initialize the number of the clusters
clusterInd = zeros(length(Xdata),2);    % initialize the array that will contain the start and end indices for each cluster
                                        % the indices are to index into Xdata and Ydata
                                        % since we don't know how large M will get, we will
                                        % initialize clusterInd to have the same rows as Xdata
clusterInd(1,1) = 1; % setting the start index for the first cluster to be 1
for i = 2:size(data,1)
    d = sqrt((Xdata(i)-Xdata(i-1))^2 + (Ydata(i)-Ydata(i-1))^2);
    if d > tol
        % start new cluster
        clusterInd(K,2) = i-1;
        K = K + 1;
        clusterInd(K,1) = i;
    end
end
clusterInd = clusterInd(clusterInd(:,1)>0,:);   % this clips off all the excess rows that only contains 0s
clusterInd(end,2) = size(data,1);               % this sets the last index to be the last element in data
%% Now that we have the clusters, let's fit the lines
lineCoeff = zeros(K,3); % initialize the vector to store the coefficients for the K line equations
lse = zeros(K,1);       % initialize the vector to store the error of the K line fits
% main line fitting code begins here
n=2;
for i = 1:K
    q = clusterInd(i,2)-clusterInd(i,1)+1;
    M = zeros(q,3);
    M(:,1) = Xdata(clusterInd(i,1):clusterInd(i,2));
    M(:,2) = Ydata(clusterInd(i,1):clusterInd(i,2));
    M(:,3) = 1;
    [U,S,V] = svd(M);
    lineCoeff(i,:) = V(:,3)';
    lse = S(3,3);
    n = n+1;
end

% the following plots each of the lines
for i = 1:K
    % computes the end points for each line segment
    xMax = max(Xdata(clusterInd(i,1):clusterInd(i,2)));
    xMin = min(Xdata(clusterInd(i,1):clusterInd(i,2)));
    if abs(lineCoeff(i,2)) < 1e-4
        xMax = -lineCoeff(i,3)/lineCoeff(i,1);
        xMin = xMax;
        yMax = max(Ydata(clusterInd(i,1):clusterInd(i,2)));
        yMin = min(Ydata(clusterInd(i,1):clusterInd(i,2)));
    else
        yMax = -lineCoeff(i,1)*xMax/lineCoeff(i,2) - lineCoeff(i,3)/lineCoeff(i,2);
        yMin = -lineCoeff(i,1)*xMin/lineCoeff(i,2) - lineCoeff(i,3)/lineCoeff(i,2);
    end
    plot([xMin xMax], [yMin yMax], 'k-')
    xMin = 0;
    xMax = 0;
    yMin = 0;
    yMax = 0;
end
hold off