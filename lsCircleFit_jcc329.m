clear all
close all
clc

% reads in the laser range finder data 
% the first column of the data denotes the angle
% the second column of the data denotes the range or distance in mm
data = load('circData.txt');

% converts the data in polar coordinates to X and Y coordinates
Xdata = data(:,2).*cos(data(:,1));
Ydata = data(:,2).*sin(data(:,1));
xm = mean(Xdata);
ym = mean(Ydata);

% plots the data
figure(1)
plot(Xdata, Ydata, 'r+')
axis([-750 1250 -1000 1000])
hold on

%% We will first separate the data set into M different clusters.  We will
%% then try to fit circles to each individual cluster

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

%% Now that we have the clusters, let's fit the circles
circCoeff = zeros(K,3); % initialize the vector to store the coefficients for the M line equations
lse = zeros(K,1);       % initialize the vector to store the error of the M line fits
% main line fitting code begins here
for i = 1:K
    q = clusterInd(i,2)-clusterInd(i,1)+1;
    u = Xdata(clusterInd(i,1):clusterInd(i,2));
    v = Ydata(clusterInd(i,1):clusterInd(i,2));
    ui = u - mean(u);
    vi = u - mean(v);
    ui2 = sum(ui.^2);
    ui3 = sum(ui.^3);
    vi2 = sum(vi.^2);
    vi3 = sum(vi.^3);
    uivi = sum(ui)*sum(vi);
    uivi2 = sum(ui.*vi.^2);
    viui2 = sum(ui.^2)*sum(vi);
    Su = sum(ui);
    Suu = sum(ui2);
    Suuu = sum(ui3);
    Sv = sum(vi);
    Svv = sum(vi2);
    Svvv = sum(vi3);
    Suv = sum(uivi);
    Suvv = sum(uivi2);
    Svuu = sum(viui2);
    uc = (.5*(Suv*(Svvv+Svuu)-(Suuu+Suvv)*Svv))/(Suu*Svv-Suv^2);
    vc = (.5*(Suvv+Svuu)-uc*Suv)/(Svv);
    alpha = uc^2+vc^2+(Suu+Svv)/q;
    xc = uc+xm;
    yc = vc+ym;
    circCoeff(i,:) = [xc;yc;alpha];
end
% the following plots each of the lines
for i = 1:K
    fprintf('I am here %d \n', i);
    theta = linspace(0,2*pi,50);
    x = circCoeff(i,1) + circCoeff(i,3)*cos(theta);
    y = circCoeff(i,2) + circCoeff(i,3)*sin(theta);
    plot(x, y, 'k-')
end
hold off