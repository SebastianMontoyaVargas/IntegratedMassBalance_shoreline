%% Read data
data_folder = '.\DUCK_TIDES\'; % Creating data list
data = dir(data_folder);
tref = datenum('01-01-1970-00:00');% Reference date
%% Maine code
% Interpolation dates
tide_time = [];
tide_level = [];

for i = 3:length(data)-1
    file_name = strcat(data_folder,data(i).name); % file name string
    time = ncread(file_name,'time'); % Reading data
    waterLevel = ncread(file_name,'waterLevel'); % Extracting water level
    time = time(~isnan(waterLevel)); % Available data dates
    waterLevel = waterLevel(~isnan(waterLevel)); % Available data
    tide_time = [tide_time; time/(3600*24)]; % Converting time to days and adding to array
    tide_level = [tide_level; waterLevel];  % Adding data to array
end
tide_time = tide_time + tref; % converting time to serial date number
tide_t = [tide_time(1):1/24:tide_time(end)]; % Interpolation date numbers
tide =  interp1(tide_time,tide_level,tide_t); % Interpolated tide
%% Saving relevant data
save('hs.mat','tide_t','tide')

