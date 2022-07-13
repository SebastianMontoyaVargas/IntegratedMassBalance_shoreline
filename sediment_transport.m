%% Read data
data_folder = '.\DUCK_DEM';
data = dir(data_folder);
load mean_waves.mat;
load Xs_DEM.mat
%% Initial parameters
gamma = 0.78;
g = 9.8;
K1 = 500;
beta = 5/12;
K2 = K1*(((g*gamma)^0.5)/2*pi())^0.2; 
%% Wave parameters
Dir = Dirm -71; % Wave direction relative to mean beach orientation
phi = deg2rad(Dir); % Wave direction in radians
% Initialize variables
dQ_dy = [];
t_f = [];
y = l_shore1;
Qy = NaN*zeros(size(Xs));
dQydy = NaN*zeros(size(Xs));
%% Data processing year
for i = 1:length(tint)
    xs = smooth(Xs(:,i)); % Smooth long-shore profile for numerical stability
    xs(isnan(Xs(:,i))) = NaN; % Eliminating false data
    if sum(isnan(xs)) < 0.2*length(xs)
        %% Long-shore slope
        for j = 1:length(xs)
            if j == 1
                dxs(j) = (xs(j+1) - xs(j))/(y(j+1) - y(j));
            elseif j == length(xs)
                dxs(j) = (xs(j) - xs(j-1))/(y(j) - y(j-1));
            else
                dxs(j) = (xs(j+1) - xs(j-1))/(y(j+1) - y(j-1));
            end
        end
        the_0 = phi(i) - atan(dxs); % wave angle relative to shoreline
        Qy(:,i) = ((K2*Hsm(i)^(12/5)*Tm(i)^(1/5)*((cos(the_0).^6).^(1/5)).*sin(the_0))'); % Sediment transport
        % Long-shore transport gradients
        for j = 1:length(y)
            if j == 1
                dQydy(j,i) = (Qy(j+1,i) - Qy(j,i))/(y(j+1) - y(j));
            elseif j == length(xs)
                dQydy(j,i) = (Qy(j,i) - Qy(j-1,i))/(y(j) - y(j-1));
            else
                dQydy(j,i) = (Qy(j+1,i) - Qy(j-1,i))/(y(j+1) - y(j-1));
            end
        end
    end
end
% Save data
save('LST.mat','tint','dQydy','Qy');