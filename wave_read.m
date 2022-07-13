%% Read data
data_folder = '.\26m_waverider';
data = dir(data_folder);

%% Data processing
WAVES_DUCK = []; % initializing array
for i = 3:length(data)
    file_name = strcat(data_folder,'\',data(i).name); % Obtaining file name
    time = ncread(file_name,'time'); % Reading time
    Hs = ncread(file_name,'waveHs'); % Reading significant wave height
    Tp = ncread(file_name,'waveTp'); % Reading peak period
    Dir = ncread(file_name,'waveMeanDirection'); % Reading wave direction
    t = []; % initializing time vector
    for j = 1:length(time)
        a = strcat('01-Jan-1970 00:00:',num2str(time(j))); % Reading time
        t = [t; datenum(a)]; % Turning time into serial time number
    end
    WAVES_DUCK = [WAVES_DUCK; t, Hs, Tp, Dir]; % Adding data to array
end
% Save data
save('WAVES_RAW.mat','WAVES_DUCK')