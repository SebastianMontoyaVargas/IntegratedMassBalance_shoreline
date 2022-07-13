%% Read data
data_folder = '.\DUCK_DEM\';
data = dir(data_folder);
tref = datenum('01-01-1970-00:00');
%% Load variables 
load('.\Xc_DEM.mat');
load('.\Xs_DEM.mat');
load('.\hs.mat')
load('.\hc_DEM.mat');
%% Interpolation coordinates
x_int = [50:.5:950]; % For DUCK DEM
Y = [-100:25:1100]; % X coordinate in Duck DEM
yvol = l_shore1;
tvol = tint;
%% Data processing
date = []; % Initialize array
for i = 4:length(data)-1
    file_name1 = strcat(data_folder,data(i-1).name);
    t = ncread(file_name1,'time');
    date = [date;  t/(3600*24)];
end
date = tref + date; 

Dc = NaN*zeros(51,length(tint)); % initialize closure depth
V1 = NaN*zeros(51,length(tint)); % initialize volume of water

%% Smoothing data
% Smoothing the influence of storms 
for j = 1:length(tint)
    Xc_smooth(:,j) = smooth(Xc(:,j),3); 
    hc_smooth(:,j) = smooth(hc(:,j),3); 
end
Xc_smooth(isnan(Xc)) = NaN; 
hc_smooth(isnan(hc)) = NaN;
%% Computing active zone width
L = Xc_smooth - Xs;
%% Computing Volumes and Dc
for i = 1:length(tint)
    t1 = date(date<=tint(i)); 
    t1 = t1(end); % Initial time for interpolation
    t2 = date(date>tint(i));
    t2 = t2(1); % Final time for interpolation
    indx1 = find(date==t1);
    indx2 = find(date==t2);
    file_name1 = strcat(data_folder,data(2+indx1(1)).name); % Data for t1
    file_name2 = strcat(data_folder,data(2+indx2(1)).name); % Data for t2
    x_shore1 = ncread(file_name1,'xFRF'); % Cross-shore position
    l_shore1 = ncread(file_name1,'yFRF'); % Long-shore position
    x_shore2 = ncread(file_name2,'xFRF');
    l_shore2 = ncread(file_name2,'yFRF');
    Z1 = ncread(file_name1,'elevation'); % Elevation at t1
    Z2 = ncread(file_name2,'elevation'); % Elevation at t2
    [X,Y,T] = ndgrid(x_shore1,l_shore1,[t1, t2]);
    Z(:,:,1) = Z1;
    Z(:,:,2) = Z2;
    Zint = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i)*ones(size(X(:,:,1)))); % Interpolated bathymetry
    hs = interp1(tide_t,smooth(tide,24*15),tint(i)); % Filtered water level at interpolation time
    Dc(:,i) = hs - hc_smooth(:,i); % Closure depth
    % Integrating water volumes
    for j = 1:length(l_shore1)
        z = interp1(x_shore1,Zint(:,j),x_int); % Interpolated profile
        xs = Xs(j,i); % Shoreline at profile
        xc = Xc_smooth(j,i); % Closure at profile
        xvol = x_int(find((x_int>xs).*(x_int<xc))); % x coordinates for integration
        zvol = z(find((x_int>xs).*(x_int<xc))); % Elevations for integration
        V1(j,i) = trapz(xvol,zvol); % Volume of water
    end
end
V = L.*Dc + V1; % Volume of sediments
beta = V./(L.*Dc); % Shape parameter
% Save data
save('vol_DEM.mat','V','beta','Dc','L')