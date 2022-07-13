clear all
load vol_DEM.mat
load Xs_DEM.mat
load mean_waves.mat
load hs.mat
load LST.mat
load Derivatives.mat
%% Deffining coordinates for south and north beach
y  = l_shore1; % Long-shore position
y_north = y(y<450); % Long-shore position north
y_south = y(y>550); % Long-shore position south
X_north = Xs(y<450,:); % Shoreline north
X_south = Xs(y>550,:); % Shoreline south

%% Attenuating effect from peak events on L and Dc
for i = 1:length(y)
    LDEM(i,:) = smooth(L(i,:),6);
    Dc_smooth(i,:) = smooth(Dc(i,:),6);
    dLdt(i,:) = smooth(dLdt(i,:),6);
    dDcdt(i,:) = smooth(dDcdt(i,:),6); 
end
Dc_smooth(isnan(Dc)) = NaN; % Removing fake data
LDEM(isnan(L)) = NaN;
%% Long-shore average and 3 std. dev boundaries
% North shoreline
XN_mean = nanmean(X_north,1); 
XN_min = nanmean(X_north,1)-3*std(X_north,1);
XN_max = nanmean(X_north,1)+3*std(X_north,1);
% South shoreline
XS_mean = nanmean(X_south,1);
XS_min = nanmean(X_south,1)-3*std(X_south,1);
XS_max = nanmean(X_south,1)+3*std(X_south,1);
% North closure depth
Dc_mean = nanmean(Dc_smooth,1);
Dc_meanN = nanmean(Dc_smooth(y<450,:),1);
Dc_minN = nanmean(Dc_smooth(y<450,:),1)-3*nanstd(Dc_smooth(y<450,:),1);
Dc_maxN = nanmean(Dc_smooth(y<450,:),1)+3*nanstd(Dc_smooth(y<450,:),1);
% South closure depth
Dc_meanS = nanmean(Dc_smooth(y>550,:),1);
Dc_minS = nanmean(Dc_smooth(y>550,:),1)-3*nanstd(Dc_smooth(y>550,:),1);
Dc_maxS = nanmean(Dc_smooth(y>550,:),1)+3*nanstd(Dc_smooth(y>550,:),1);
Dc_min = min(Dc_smooth); % Auxiliar variable
Dc_max = max(Dc_smooth);

% Active zone width
Lmean = nanmean(LDEM,1);
Lmin = min(LDEM);
Lmax = max(LDEM);
LmeanN = nanmean(LDEM(y<450,:),1);
LmeanS = nanmean(LDEM(y>550,:),1);
LminN = nanmean(LDEM(y<450,:),1)-3*nanstd(LDEM(y<450,:),1);
LminS = nanmean(LDEM(y>550,:),1)-3*nanstd(LDEM(y>550,:),1);
LmaxN = nanmean(LDEM(y<450,:),1)+3*nanstd(LDEM(y<450,:),1);
LmaxS = nanmean(LDEM(y>550,:),1)+3*nanstd(LDEM(y>550,:),1);

% Shape parameter
betamean = nanmean(beta,1);
betameanN = nanmean(beta(y<450,:),1);
betameanS = nanmean(beta(y>550,:),1);
betamin = min(beta);
betamax = max(beta);
betaminN = nanmean(beta(y<450,:),1)-3*nanstd(beta(y<450,:),1);
betaminS = nanmean(beta(y>550,:),1)-3*nanstd(beta(y>550,:),1);
betamaxN = nanmean(beta(y<450,:),1)+3*nanstd(beta(y<450,:),1);
betamaxS = nanmean(beta(y>550,:),1)+3*nanstd(beta(y>550,:),1);

%% North and south shoreline changes
XdotN = dXdt(y<450,:); % North
XdotS = dXdt(y>550,:); % South

%% Shoreline change model terms
Xdot1 = -beta.*dLdt;
Xdot2 = (1-beta).*LDEM./Dc_smooth.*dDcdt;
Xdot3 = -LDEM.*dbetadt;
Xdot4 = -LDEM./Dc_smooth.*(ones(51,1)*dhsdt);
Xdot5 = -1./Dc_smooth.*dQydy;
Xdot_model = Xdot1+Xdot2+Xdot3+Xdot4+Xdot5;

%% Cross-correlation
% initialize variables
R1_S = NaN*[1:48]'; % 1st term correlation south
R2_S = NaN*[1:48]'; % 2nd term correlation south
R3_S = NaN*[1:48]'; % 3rd term correlation south
R4_S = NaN*[1:48]'; % 4th term correlation south
R5_S = NaN*[1:48]'; % 5th term correlation south
R1_N = NaN*[1:48]'; % 1st term correlation north
R2_N = NaN*[1:48]'; % 2nd term correlation north
R3_N = NaN*[1:48]'; % 3rd term correlation north
R4_N = NaN*[1:48]'; % 4th term correlation north
R5_N = NaN*[1:48]'; % 5th term correlation north
R12_N = NaN*[1:48]'; % Term 1+2 north
R12_S = NaN*[1:48]'; % term 1 + 2 south
R123_N = NaN*[1:48]'; % term 1+2+3 north
R123_S = NaN*[1:48]'; % term 1+2+3 south
% Computing correlations
for i = 1:48
    % Initialize auxiliar variables
    Xdot_filt1 = zeros(size(Xdot1));
    Xdot_filt2 = zeros(size(Xdot2));
    Xdot_filt3 = zeros(size(Xdot3));
    Xdot_filt4 = zeros(size(Xdot4));
    Xdot_filt5 = zeros(size(Xdot5));
    % Processing time filtering for each cross-section
    for j = 1:length(y)
        Xdot_filt1(j,:) = smooth(Xdot1(j,:),i);
        Xdot_filt2(j,:) = smooth(Xdot2(j,:),i);
        Xdot_filt3(j,:) = smooth(Xdot3(j,:),i);
        Xdot_filt4(j,:) = smooth(Xdot4(j,:),i);
        Xdot_filt5(j,:) = smooth(Xdot5(j,:),i);
    end
    % Compute correlations neglecting blanck data
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt1(y<450,10:70),1),'Rows','Complete');
    R1_N(i) = A(1,2); 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt1(y>550,10:70),1),'Rows','Complete');
    R1_S(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt2(y<450,10:70),1),'Rows','Complete');
    R2_N(i) = [ A(1,2)]; 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt2(y>550,10:70),1),'Rows','Complete');
    R2_S(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt3(y<450,10:70),1),'Rows','Complete');
    R3_N(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt3(y>550,10:70),1),'Rows','Complete');
    R3_S(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt1(y<450,10:70)+Xdot_filt2(y<450,10:70),1),'Rows','Complete');
    R12_N(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt1(y>550,10:70)+Xdot_filt2(y>550,10:70),1),'Rows','Complete');
    R12_S(i) = [A(1,2)];  
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt1(y<450,10:70)+Xdot_filt2(y<450,10:70)+Xdot_filt3(y<450,10:70),1),'Rows','Complete');
    R123_N(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt1(y>550,10:70)+Xdot_filt2(y>550,10:70)+Xdot_filt3(y>550,10:70),1),'Rows','Complete');
    R123_S(i) = [ A(1,2)]; 
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt4(y<450,10:70),1),'Rows','Complete');
    R4_N(i) = [A(1,2)]; 
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt4(y>550,10:70),1),'Rows','Complete');
    R4_S(i) = A(1,2);
    A = corrcoef(nanmean(XdotN(:,10:70),1),nanmean(Xdot_filt5(y<450,10:70),1),'Rows','Complete');
    R5_N(i) = A(1,2);
    A = corrcoef(nanmean(XdotS(:,10:70),1),nanmean(Xdot_filt5(y>550,10:70),1),'Rows','Complete');
    R5_S(i) = A(1,2);
end
win_size = [1:48]*15;
save('Correlations.mat','R1_N', 'R2_N','R3_N','R4_N','R5_N','R12_N',...
    'R123_N','R1_S', 'R2_S','R3_S','R4_S','R5_S','R12_S','R123_S','win_size')