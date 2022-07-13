clear all
%% Load data
load vol_DEM.mat
load Xs_DEM.mat
load mean_waves.mat
load hs.mat
%% Main code
% Filter sea level
hs = interp1(tide_t,smooth(tide,24*15),tint);
% Initialize variables
dLdt = zeros(size(L));
dXdt = zeros(size(Xs));
dDcdt = zeros(size(Dc));
dhsdt = zeros(size(tint));
dbetadt = zeros(size(beta));
% Compute initial and last finite differences 
dLdt(:,1) = (L(:,2) - L(:,1))./(tint(2) - tint(1)); % Forward difference
dLdt(:,end) = (L(:,end) - L(:,end-1))./(tint(end) - tint(end-1)); % Backward difference
dXdt(:,1) = (Xs(:,2) - Xs(:,1))./(tint(2) - tint(1));
dXdt(:,end) = (Xs(:,end) - Xs(:,end-1))./(tint(end) - tint(end-1));
dDcdt(:,1) = (Dc(:,2) - Dc(:,1))./(tint(2) - tint(1));
dDcdt(:,end) = (Dc(:,end) - Dc(:,end-1))./(tint(end) - tint(end-1));
dbetadt(:,1) = (beta(:,2) - beta(:,1))./(tint(2) - tint(1));
dbetadt(:,end) = (beta(:,end) - beta(:,end-1))./(tint(end) - tint(end-1));
dhsdt(:,1) = (hs(2) - hs(1))./(tint(2) - tint(1));
dhsdt(:,end) = (hs(end) - hs(end-1))./(tint(end) - tint(end-1));

% Computing central difference for intermediate times
for i = 2:length(tint)-1
    dLdt(:,i) = (L(:,i+1) - L(:,i-1))./(tint(i+1) - tint(i-1));
    dXdt(:,i) = (Xs(:,i+1) - Xs(:,i-1))./(tint(i+1) - tint(i-1));
    dDcdt(:,i) = (Dc(:,i+1) - Dc(:,i-1))./(tint(i+1) - tint(i-1));
    dbetadt(:,i) = (beta(:,i+1) - beta(:,i-1))./(tint(i+1) - tint(i-1));
    dhsdt(i) = (hs(i+1) - hs(i-1))./(tint(i+1) - tint(i-1));
end
%% Saving data
save('Derivatives.mat','dXdt','dLdt','dDcdt','dbetadt','dhsdt')