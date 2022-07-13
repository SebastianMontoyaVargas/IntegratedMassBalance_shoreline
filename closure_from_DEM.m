%% Read data
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

Xc = NaN*zeros(51,length(tint)); % Initializing closure line location array
hc = NaN*zeros(51,length(tint)); % Initializing closure line location array
for i = 2:length(tint)
        t1 = date(date<=tint(i)); % Available dates before current interpolation time
        t12 = t1(end); % Date of measurement previous to the inteprolation date
        t2 = date(date>tint(i)); % Available dates after current interpolation time
        t22 = t2(1); % Date of measurement after the interpolation date
        
        t1 = date(date<=tint(i-1)); % Available dates before current interpolation time
        t11 = t1(end); % Date of measurement previous to the inteprolation date
        t2 = date(date>tint(i-1)); % Available dates after current interpolation time
        t21 = t2(1); % Date of measurement after the interpolation date
        
        if t21 == t11
            indx11 = find(date==t11); % Corresponding index of t11
            indx21 = find(date==t11); % Corresponding index of t21
            indx12 = find(date==t12); % Corresponding index of t12
            indx22 = find(date==t12); % Corresponding index of t22

            file_name11 = strcat(data_folder,data(2+indx11(1)).name); % data file for t11
            file_name21 = strcat(data_folder,data(2+indx21(1)).name); % data file for t21
            file_name12 = strcat(data_folder,data(2+indx12(1)).name); % data file for t12
            file_name22 = strcat(data_folder,data(2+indx22(1)).name); % data file for t22

            x_shore1 = ncread(file_name12,'xFRF'); % cross-shore coordinate for t11
            l_shore1 = ncread(file_name12,'yFRF'); % long-shore coordinate for t11
            Z1 = ncread(file_name12,'elevation'); % DEM for t11
            Z2 = ncread(file_name22,'elevation'); % DEM for t21
            
            [X,Y,T] = ndgrid(x_shore1,l_shore1,[t11, t21]); % 4D Interpolation grid 
            Z(:,:,1) = Z1; 
            Z(:,:,2) = Z2;
            Zint1 = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i-1)*ones(size(X(:,:,1)))); % Interpolated DEM at time i
            Zint2 = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i)*ones(size(X(:,:,1)))); % Interpolated DEM at time i
            
        else
            indx11 = find(date==t11); % Corresponding index of t11
            indx21 = find(date==t21); % Corresponding index of t21
            indx12 = find(date==t12); % Corresponding index of t12
            indx22 = find(date==t22); % Corresponding index of t22

            file_name11 = strcat(data_folder,data(2+indx11(1)).name); % data file for t11
            file_name21 = strcat(data_folder,data(2+indx21(1)).name); % data file for t21
            file_name12 = strcat(data_folder,data(2+indx12(1)).name); % data file for t12
            file_name22 = strcat(data_folder,data(2+indx22(1)).name); % data file for t22

            x_shore1 = ncread(file_name12,'xFRF'); % cross-shore coordinate for t11
            l_shore1 = ncread(file_name12,'yFRF'); % long-shore coordinate for t11
            Z1 = ncread(file_name12,'elevation'); % DEM for t11
            Z2 = ncread(file_name22,'elevation'); % DEM for t21
            [X,Y,T] = ndgrid(x_shore1,l_shore1,[t12, t22]); % 4D Interpolation grid 
            Z(:,:,1) = Z1; 
            Z(:,:,2) = Z2;
            Zint2 = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i)*ones(size(X(:,:,1)))); % Interpolated DEM at time i


            x_shore1 = ncread(file_name11,'xFRF'); % cross-shore coordinate for t11
            l_shore1 = ncread(file_name11,'yFRF'); % long-shore coordinate for t11
            Z1 = ncread(file_name11,'elevation'); % DEM for t11
            Z2 = ncread(file_name21,'elevation'); % DEM for t21
            [X,Y,T] = ndgrid(x_shore1,l_shore1,[t11, t21]); % 4D Interpolation grid 
            Z(:,:,1) = Z1; 
            Z(:,:,2) = Z2;
            Zint1 = interpn(X,Y,T,Z,X(:,:,1),Y(:,:,1),tint(i-1)*ones(size(X(:,:,1)))); % Interpolated DEM at time i-1
        end
        
        for j = 1:length(l_shore1)
            zi = interp1(x_shore1,Zint1(:,j),x_int); % Profile at coordinate y and initial time
            zf = interp1(x_shore1,Zint2(:,j),x_int); % Profile at coordinate y and final time
            dz = abs(zf-zi); % change in elevation between time step
            dz_aux = smooth(abs(dz(200:1400)-0.1),100);
            [M I] = min(dz_aux(100:end-100)); % Finding most offshore value under 0.05 m
            if sum(isnan(dz)) > length(dz)/2
                Xc(j,i) = NaN; % Insuficient data
            else
                Xc(j,i) = x_int(300+I); % Closure at the fartest profile change under 0.05 m
            end
            hc(j,i) = zf(300+I);
        end
end
% Removing outliers
hc(abs(normalize(Xc))>1) = NaN;
Xc(abs(normalize(Xc))>1) = NaN;
% Saving data
save('Xc_DEM.mat','Xc','tint','l_shore1')
save('hc_DEM.mat','hc','tint','l_shore1')