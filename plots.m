clear all
addpath('.\DATA_PROCESSING\')
load vol_DEM.mat
load Xs_DEM.mat
load mean_waves.mat
load hs.mat
load LST.mat
load WAVES_RAW.mat
load Derivatives.mat
load Correlations.mat
%% Deffining coordinates for south and north beach
y  = l_shore1; % Long-shore position
level = interp1(tide_t,smooth(tide,24),tint);
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

%% Predicted shoreline change
xdot_rel = zeros(size(Xdot1));
Xdot_filtN = zeros(size(Xdot1));
Xdot_filtS = zeros(size(Xdot1));
xdot123N = zeros(size(Xdot1));
xdot123S = zeros(size(Xdot1));
xdot5filt = zeros(size(Xdot1));
xdot4filt = zeros(size(Xdot1));

for i = 1:length(y)
    Xdot_filtN(i,:) = smooth(Xdot1(i,:)+Xdot2(i,:)+Xdot3(i,:),28) + smooth(Xdot4(i,:),4) + smooth(Xdot5(i,:),3);
    Xdot_filtS(i,:) = smooth(Xdot1(i,:)+Xdot2(i,:)+Xdot3(i,:),28) + smooth(Xdot4(i,:),4) + smooth(Xdot5(i,:),3);
    xdot123N(i,:) = (smooth(Xdot1(i,:)+Xdot2(i,:)+Xdot3(i,:),28))';
    xdot123S(i,:) = smooth(Xdot1(i,:)+Xdot2(i,:)+Xdot3(i,:),1)';
    xdot5filt(i,:) = smooth(Xdot5(i,:),3);
    xdot4filt(i,:) = smooth(Xdot4(i,:),4);
end

%% Filtered model
XdotN_filt = nanmean(Xdot_filtN(y<450,:),1);
XdotS_filt = nanmean(Xdot_filtS(y>550,:),1);
%% Maximum and minimum values for filtered series
XdotN_min = min(Xdot_filtN(y<450,:));
XdotN_max = max(Xdot_filtS(y<450,:));
XdotS_min = min(Xdot_filtN(y>550,:));
XdotS_max = max(Xdot_filtS(y>550,:));

%% Long-shore average for filtered series
Xdot_fNmean = nanmean(XdotN_filt,1);
Xdot_fSmean = nanmean(XdotS_filt,1);
Xdot_mNmean = nanmean(XdotN,1);
Xdot_4mean = nanmean(xdot4filt,1);
Xdot_mSmean = nanmean(XdotS,1);

%% Yearly correlations and wave average
% Initialize variables
RN_time = zeros(length(tint)-24,1);
RS_time = zeros(length(tint)-24,1);
Hsy = zeros(length(tint)-24,1);
Ty = zeros(length(tint)-24,1);;
Diry = zeros(length(tint)-24,1);;
% Computing correlations
for i = 1:length(tint)-24
    % North beach
    A = corrcoef(Xdot_mNmean(i:i+24),Xdot_fNmean(i:i+24),'Rows','Complete');
    RN_time(i) = A(1,2); 
    % South beach
    A = corrcoef(Xdot_mSmean(i:i+24),Xdot_fSmean(i:i+24),'Rows','Complete');
    RS_time(i) = A(1,2); 
    % Average waves
    Hsy(i) = nanmean(Hsm(i:i+24));
    Ty(i) = nanmean(Tm(i:i+24));
    Diry(i) = nanmean(Dirm(i:i+24));
end
%% Plots
figure(7)
subplot(2,2,1)
hold on;
grid on;
plot(tint(1:192),RN_time,'r')
plot(tint(1:192),RS_time,'b')
xlabel('Year')
ylabel('Corr. Coef')
datetick('x','yyyy','keeplimits')

subplot(2,2,3)
hold on;
grid on;
plot(tint(1:192),Hsy./(9.81.*Ty.^2),'k')
ylabel('Yearly Avg. H_s/(gT^2)')
yyaxis right
plot(tint(1:192),Diry-71)
xlabel('Year')
ylabel('Yearly Avg. \theta_0 (°)')
datetick('x','yyyy','keeplimits')

subplot(2,2,2)
hold on;
grid on;
title('Corr. Coef. North')
scatter((Hsy./(9.81.*Ty.^2)),Diry-71,25,RN_time,'filled')
colorbar
xlabel('Yearly Avg. H_s/(gT^2)')
ylabel('Yearly Avg. \theta_0 (°)')
subplot(2,2,4)
hold on;
grid on;
title('Corr. Coef. South')
scatter((Hsy./(9.81.*Ty.^2)),Diry-71,25,RS_time,'filled')
colorbar
xlabel('Yearly Avg. H_s/(gT^2)')
ylabel('Yearly Avg. \theta_0 (°)')


%%----------------
figure(6)
subplot(2,1,1)
hold on;
grid on;
plot(tint,(XdotN(1,:)),'.','Color',0.8*[1 1 1])
plot(tint,nanmean(XdotN),'r','LineWidth',1)
plot(tint,nanmean(Xdot_filtN),'k','LineWidth',1)
plot(tint,(XdotN),'.','Color',0.8*[1 1 1])
plot(tint,nanmean(XdotN),'r','LineWidth',1)
plot(tint,nanmean(Xdot_filtN),'k','LineWidth',1)
datetick('x','yyyy','keeplimits')
ylabel('North \partialX_s/\partialt (m/day)')
ylim([-2 2])
legend('DEM Meas.','Long-shore avg. Meas.','Long-shore avg. Proposed Model')

subplot(2,1,2)
hold on;
grid on;
plot(tint,(XdotS),'.','Color',0.8*[1 1 1])
plot(tint,nanmean(XdotS),'r','LineWidth',1)
plot(tint,nanmean(Xdot_filtS),'k','LineWidth',1)
datetick('x','yyyy','keeplimits')
ylabel('South \partialX_s/\partialt (m/day)')
xlabel('Year')
ylim([-2 2])

% ----------------------------
figure(5)
hold on
grid on
plot(LDEM(:),Dc_smooth(:),'.','Color',[0.7 0.7 0.7])
plot([0:600],0.195*[0:600].^0.56,'k','LineWidth',1.3)
plot([0:600],0.108*[0:600].^(2/3),'r','LineWidth',1.3)
xlabel('Active zone width, L (m)')
ylabel('Closure depth, D_c (m)')
xlim([0 600])
legend('Measurement','D_c = 0.195*L^{0.56}','D_c = 0.108*L^{2/3}','Location','southeast')
% ----------------------------
figure(4)
subplot(4,1,1)
hold on
plot(tide_t,tide,'.','Color',[0.7 0.7 0.7])
plot(tint,level,'k','LineWidth',1.3)
plot(tint,smooth(level,4),'r','LineWidth',1.1)
datetick('x','yyyy','keeplimits')
xlim([tint(1) tint(end)])
ylabel('h_s (m)')
legend('Hourly','15-day average','60-day average','Orientation','horizontal' )
subplot(4,1,2)
hold on
plot(WAVES_DUCK(:,1),WAVES_DUCK(:,2),'.','Color',[0.7 0.7 0.7])
plot(tint,Hsm,'k','LineWidth',1.3)
plot(tint,smooth(Hsm,24),'r','LineWidth',1.1)
xlim([tint(1) tint(end)])
datetick('x','yyyy','keeplimits')
ylabel('H_s (m)')
legend('Hourly','15-day average','360-day average','Orientation','horizontal')
subplot(4,1,3)
hold on
plot(WAVES_DUCK(:,1),WAVES_DUCK(:,3),'.','Color',[0.7 0.7 0.7])
plot(tint,Tm,'k','LineWidth',1.3)
plot(tint,smooth(Tm,24),'r','LineWidth',1.1)
xlim([tint(1) tint(end)])
datetick('x','yyyy','keeplimits')
ylabel('T_p (m)')
legend('Hourly','15-day average','360-day average','Orientation','horizontal')
subplot(4,1,4)
hold on
plot(WAVES_DUCK(:,1),WAVES_DUCK(:,4),'.','Color',[0.7 0.7 0.7])
plot(tint,Dirm,'k','LineWidth',1.3)
plot(tint,smooth(Dirm,24),'r','LineWidth',1.1)
xlim([tint(1) tint(end)])
ylim([0 360])
datetick('x','yyyy','keeplimits')
ylabel('\phi_0 (º)')
xlabel('Year')
yticks([0 90 180 270])
yyaxis right
ylim([0 360])
plot(tide_t,tide,'LineStyle','none')
yticks([0 90 180 270])
yticklabels({'N','E','S','W'})
legend('Hourly','15-day average','360-day average','Orientation','horizontal')
% -------------------------
figure(3)
subplot(6,1,1)
hold on
plot(win_size,R1_N,'k')
plot(win_size,R1_S,'b')
title('-\beta*\partialL/\partialt')
legend('North','South')
ylabel('Corr.Coef')
% ylim([-.1 0.4])

subplot(6,1,2)
hold on
plot(win_size,R2_N,'k')
plot(win_size,R2_S,'b')
title('(1-\beta)*L/D_c*\partialD_c/\partialt')
% legend('North','South')
ylabel('Corr.Coef')
% ylim([-.40 0.2])


subplot(6,1,3)
hold on
plot(win_size,R3_N,'k')
plot(win_size,R3_S,'b')
title('-L*\partial\beta/\partialt')
% legend('North','South')
ylabel('Corr.Coef')
% ylim([-.2 0.2])



subplot(6,1,4)
hold on
plot(win_size,R123_N,'b')
plot(win_size,R123_S,'r')
title('-\beta*(\partialL/\partialt)_{CS}+(1-\beta)*L/D_c*\partialD_c/\partialt-L*\partial\beta/\partialt')
% legend('North','South')
ylabel('Corr.Coef')
% xlabel('Window size (days)')
% ylim([-.1 0.45])


subplot(6,1,5)
hold on
plot(win_size,R4_N,'k')
plot(win_size,R4_S,'b')
title('-L/D_c*\partialh_s/\partialt')
% legend('North','South')
ylabel('Corr.Coef')
% xlabel('Window size (days)')
% ylim([-.35 0.65])
yticks([-.3 0 .3 .6])

subplot(6,1,6)
hold on
plot(win_size,R5_N,'k')
plot(win_size,R5_S,'b')
title('-1/D_c*\partialQ_y/\partialy')
% legend('North','South')
ylabel('Corr.Coef')
xlabel('Window size (days)')


%%------------------------------------------------------
figure(2)
subplot(5,1,1)
hold on
plot(tint,nanmean(Xdot1(y<450,:),1),'k','LineWidth',1)
plot(tint,nanmean(Xdot1(y<450,:),1)-nanstd(Xdot1(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot1(y<450,:),1)+nanstd(Xdot1(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot1(y>550,:),1),'b','LineWidth',1)
plot(tint,nanmean(Xdot1(y>550,:),1)-nanstd(Xdot1(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
plot(tint,nanmean(Xdot1(y>550,:),1)+nanstd(Xdot1(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
ylim([-5 5])
ylabel(['-\beta*(\partialL/\partialt)' newline '(m/day)'])
% legend('L.S. avg.','-\beta*(\partialL/\partialt)_{CS}')
datetick('x','yyyy','keeplimits')

subplot(5,1,2)
hold on
plot(tint,nanmean(Xdot2(y<450,:),1),'k','LineWidth',1)
plot(tint,nanmean(Xdot2(y<450,:),1)-nanstd(Xdot2(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot2(y<450,:),1)+nanstd(Xdot2(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot2(y>550,:),1),'b','LineWidth',1)
plot(tint,nanmean(Xdot2(y>550,:),1)-nanstd(Xdot2(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
plot(tint,nanmean(Xdot2(y>550,:),1)+nanstd(Xdot2(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
ylim([-5 5])
ylabel(['(1-\beta)*L/D_c*(\partialD_c/\partialt)' newline ' (m/day)'])
% legend('L.S. avg.','-\beta*(\partialL/\partialt)_{CS}')
datetick('x','yyyy','keeplimits')


subplot(5,1,3)
hold on
plot(tint,nanmean(Xdot3(y<450,:),1),'k','LineWidth',1)
plot(tint,nanmean(Xdot3(y<450,:),1)-nanstd(Xdot3(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot3(y<450,:),1)+nanstd(Xdot3(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot3(y>550,:),1),'b','LineWidth',1)
plot(tint,nanmean(Xdot3(y>550,:),1)-nanstd(Xdot3(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
plot(tint,nanmean(Xdot3(y>550,:),1)+nanstd(Xdot3(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
% ylim([-5 5])
ylabel(['-L*(\partial\beta/\partialt)' newline ' (m/day)'])
% legend('L.S. avg.','-\beta*(\partialL/\partialt)_{CS}')
datetick('x','yyyy','keeplimits')

subplot(5,1,4)
hold on
plot(tint,nanmean(Xdot4(y<450,:),1),'k','LineWidth',1)
plot(tint,nanmean(Xdot4(y<450,:),1)-nanstd(Xdot4(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot4(y<450,:),1)+nanstd(Xdot4(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot4(y>550,:),1),'b','LineWidth',1)
plot(tint,nanmean(Xdot4(y>550,:),1)-nanstd(Xdot4(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
plot(tint,nanmean(Xdot4(y>550,:),1)+nanstd(Xdot4(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
% ylim([-5 5])
ylabel(['-L/D_c*(\partialh_s/\partialt)' newline ' (m/day)'])
% legend('L.S. avg.','-\beta*(\partialL/\partialt)_{CS}')
datetick('x','yyyy','keeplimits')


subplot(5,1,5)
hold on
plot(tint,nanmean(Xdot5(y<450,:),1),'k','LineWidth',1)
plot(tint,nanmean(Xdot5(y<450,:),1)-nanstd(Xdot5(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot5(y<450,:),1)+nanstd(Xdot5(y<450,:),1),'--','Color',[0.5 0.5 0.5])
plot(tint,nanmean(Xdot5(y>550,:),1),'b','LineWidth',1)
plot(tint,nanmean(Xdot5(y>550,:),1)-nanstd(Xdot5(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
plot(tint,nanmean(Xdot5(y>550,:),1)+nanstd(Xdot5(y>550,:),1),'--','Color',0.5*[0 0.5 0.7])
% ylim([-5 5])
ylabel(['-1/D_c*(\partialQ_y/\partialy) ' newline '(m/day)'])
% legend('L.S. avg.','-\beta*(\partialL/\partialt)_{CS}')
datetick('x','yyyy','keeplimits')

xlabel('Year')

%%-------------------------------
figure(1)
subplot(4,1,1)
plot(tint,XN_mean,'k','LineWidth',1)
hold on;
plot(tint,XN_min,'--k')
plot(tint,XN_max,'--k')
plot(tint,XS_mean,'b','LineWidth',1)
plot(tint,XS_min,'--b')
plot(tint,XS_max,'--b')
datetick('x','yyyy','keeplimits')
ylabel('x_s (m)')
% legend('North','South','Orientation','horizontal' )
xlim([tint(1) tint(end)])

subplot(4,1,2)
plot(tint,Dc_meanN,'k','LineWidth',1)
hold on;
plot(tint,Dc_meanS,'b','LineWidth',1)
plot(tint,Dc_minN,'--k')
plot(tint,Dc_maxN,'--k')
plot(tint,Dc_minS,'--b')
plot(tint,Dc_maxS,'--b')
datetick('x','yyyy','keeplimits')
ylabel('D_c (m)')
% legend('Mean','Max./Min.','Orientation','horizontal' )
xlim([tint(1) tint(end)])

subplot(4,1,3)
plot(tint,LmeanN,'k')
hold on;
plot(tint,LmeanS,'b')
plot(tint,LminN,'--k')
plot(tint,LmaxN,'--k')
plot(tint,LminS,'--b')
plot(tint,LmaxS,'--b')
datetick('x','yyyy','keeplimits')
% xlabel('Year')
ylabel('L (m)')
% legend('Mean','Max. and Min.','Orientation','horizontal','Orientation','horizontal'  )
xlim([tint(1) tint(end)])

subplot(4,1,4)
plot(tint,betameanN,'k')
hold on;
plot(tint,betameanS,'b')
plot(tint,betaminN,'--k')
plot(tint,betamaxN,'--k')
plot(tint,betaminS,'--b')
plot(tint,betamaxS,'--b')
datetick('x','yyyy','keeplimits')
xlabel('Year')
ylabel('\beta ')
% legend('Mean','Max. and Min.','Orientation','horizontal','Orientation','horizontal'  )
xlim([tint(1) tint(end)])