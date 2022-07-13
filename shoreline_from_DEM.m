%% Read data
load '.\hs.mat';
data_folder = '.\DUCK_DEM\';
data = dir(data_folder); % Creating data list
tref = datenum('01-01-1970-00:00'); % Reference date
%% Maine code

% Creating date-time array
date = []; % Initializing date-time array
for i = 4:length(data)-1
    file_name1 = strcat(data_folder,data(i-1).name); % creating file name
    t = ncread(file_name1,'time'); % reading time in seconds from first reading
    date = [date;  t/(3600*24)]; % converting time into days and adding to array
end
date = tref + date; % Set serial number dates to the reference date number

% Interpolation
x_int = [50:.5:950]; % x coordinate for interpolation 
tint = [date(276):15:date(end)]; % time for interpolation (15 days intervals)

Xs = NaN*zeros(51,length(tint)); % Initializing shoreline location array

for i = 1:length(tint)
        t1 = date(date<=tint(i)); % Available dates before current interpolation time
        t1 = t1(end); % Date of measurement previous to the inteprolation date
        t2 = date(date>tint(i)); % Available dates after current interpolation time
        t2 = t2(1); % Date of measurement after the interpolation date
        indx1 = find(date==t1); % Corresponding index of t1
        indx2 = find(date==t2); % Corresponding index of t2
        file_name1 = strcat(data_folder,data(2+indx1(1)).name); % data file for t1
        file_name2 = strcat(data_folder,data(2+indx2(1)).name); % data file for t2
        x_shore1 = ncread(file_name1,'xFRF'); % cross-shore coordinate for t1
        l_shore1 = ncread(file_name1,'yFRF'); % long-shore coordinate for t1
        Z1 = ncread(file_name1,'elevation'); % DEM for t1
        Z2 = ncread(file_name2,'elevation'); % DEM for t2
        [X,Y,T] = ndgrid(x_shore1,l_shore1,[t1, t2]); % 4D Interpolation grid 
        Z(:,:,1) = Z1; 
        Z(:,:,2) = Z2;
        Zint = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i)*ones(size(X(:,:,1)))); % Interpolated DEM
        water_level = interp1(tide_t,smooth(tide,24*15),tint(i)); % Filtered water level at interpolation time
        for j = 1:length(l_shore1)
            z = interp1(x_shore1,Zint(:,j),x_int); % Profile at coordinate y
            dz = abs(z-water_level); % Elevation relative to water level
            [M I] = min(dz); % Finding nearest zero crossing
            if sum(isnan(dz)) > length(dz)/2
                Xs(j,i) = NaN; % Blank data
            else
                Xs(j,i) = x_int(I); % shoreline position as intersection between water level and beach profile
            end
        end
end
save('Xs_DEM.mat','Xs','tint','l_shore1')