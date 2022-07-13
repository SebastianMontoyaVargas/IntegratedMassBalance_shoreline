%% Load data
load 'WAVES_RAW.mat';
load('Xs_DEM','tint')

Hst = WAVES_DUCK(:,1);
Hs = WAVES_DUCK(:,2);
T =  WAVES_DUCK(:,3);
Dir =  WAVES_DUCK(:,4);
%% Initialize variables
Hsm = []; % Mean Hs
Tm = []; % Mean Tp
Dirm = []; % Mean Dir

% Main code
for i = 1:length(tint)
    if i <length(tint)
        % Set of data whitin time interval
        wave = Hs(Hst>tint(i)&Hst<tint(i+1)); 
        period = T(Hst>tint(i)&Hst<tint(i+1));
        direction = Dir(Hst>tint(i)&Hst<tint(i+1));
        % Average values
        Hsmean = nanmean(wave);
        Tmean = nanmean(period);
        Dirmean = nanmean(direction);
    else
        % Set of data for last interpolation time
        auxt = Hst>tint(i);
        wave = Hs(auxt(1:24*15));
        period = T(auxt(1:24*15));
        direction = Dir(auxt(1:24*15));
        % Average
        Hsmean = nanmean(wave);
        Tmean = nanmean(period);
        Dirmean = nanmean(direction);
    end
    %% Adding data to arrays
    Hsm = [Hsm; Hsmean];
    Tm = [Tm; Tmean];
    Dirm = [Dirm; Dirmean];
end
% Saving data
save('mean_waves.mat','tint','Hsm','Tm','Dirm')