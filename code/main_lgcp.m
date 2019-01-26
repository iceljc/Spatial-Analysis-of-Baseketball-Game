clc;clear;

%%%%%%%%  shot plot  %%%%%%%%%%%%
%%
f = xlsread('nba_savantHarden.csv');

label = f(:, 8);
x = f(:, 13);
y = f(:, 14);
shot_distance = f(:, 11);


index_long_dist = find(shot_distance >= 35);
label(index_long_dist) = [];
x(index_long_dist) = [];
y(index_long_dist) = [];

index_shotmade = find(label == 1);
index_shotmiss = find(label == 0);

x_made = x(index_shotmade);
y_made = y(index_shotmade);

x_miss = x(index_shotmiss);
y_miss = y(index_shotmiss);

figure
plot(x_made, y_made, 'ob', 'MarkerSize', 10, 'LineWidth', 3); hold on;
plot(x_miss, y_miss, 'xr', 'MarkerSize', 10, 'LineWidth', 3);
 
% N = 50;
% x = linspace(-250, 250, N);
% y = linspace(-50, 350, N);
% plot(x,meshgrid(y,x),'k'); hold on
% plot(meshgrid(x,y),y,'k'); 
% axis off;
court_plot();
pause(2)
close all;

%% count matrix and plot

N = 50;
xd = linspace(-250, 250, N);
yd = linspace(-50, 350, N);

xc = linspace(1, N-1, N-1);
yc = linspace(1, N-1, N-1);

% xc = linspace(1 - N/2, N/2 - 1, N-1);
% yc = linspace(1 - N/2, N/2 - 1, N-1);

xq = x;
yq = y;

count = zeros(N-1, N-1);

for j = 1 : N-1
    for i = 1 : N-1
        xv = [xd(i), xd(i+1), xd(i+1), xd(i), xd(i)];
        yv = [yd(j), yd(j), yd(j+1), yd(j+1), yd(j)];
%         center_x(i,j) = (xd(i) + xd(i+1)) / 2;
%         center_y(i,j) = (yd(i) + yd(i+1)) / 2;        
        center_x(i,j) = xc(i);
        center_y(i,j) = yc(i); 
        in = inpolygon(xq, yq, xv, yv);
        count(i,j) = sum(in(:));

    end
end

cx = reshape(center_x, (N-1) * (N-1), 1);
cy = reshape(center_y', 1, (N-1) * (N-1));
% cy = fliplr(cy);

center(1, :) = cx(:, 1);
center(2, :) = cy(1, :);

count = rot90(count);



figure
imagesc(count); colorbar('FontSize', 15); %colormap(flipud(autumn));
set(gca,'xtick',[], 'ytick', []);
set(gcf,'position',[200,50,800,600])
pause(2)
close all;

%% LGCP surface


% sigma = std(count(:));
variance = 1;

Xv = reshape(count', 1, (N-1)*(N-1));

save main Xv

z0 = 0; % bias
phi = 4; % length scale

K = zeros((N-1)*(N-1)); % covariance matrix

for i = 1 : (N-1)*(N-1)
    for j = 1 : (N-1)*(N-1)

        K(i,j) = variance * exp(-0.5 * (center(:,i) - center(:,j))'*(center(:,i) - center(:,j)) / phi^2);
    end
end

 
K_mod = K + 0.00000001*eye(size(K)); % modified covariance matrix


z_init = 0.0 * ones(1, (N-1)*(N-1));


prior = chol(K_mod, 'lower'); % prior

[zn, cur_log_like] = elliptical_slice(z_init, prior, @log_likelihood);


k = 1; % iteration number

while 1
    
    z_init = zn;
    [zn, cur_log_like] = elliptical_slice(z_init, prior, @log_likelihood);
    k = k + 1;
    
    convcheck = norm(zn - z_init, 'fro') / norm(zn, 'fro');
    
    if convcheck < 1e-5 %|| k == 200
        break;
    end
    
end
    



lamda = exp(zn + z0);

dA = (xv(2)-xv(1))*(yv(3)-yv(1));

lamda = reshape(lamda', N-1, N-1);


figure
imagesc(lamda'); colorbar; %colormap(flipud(autumn));
set(gca,'xtick',[], 'ytick', []);
set(gcf,'position',[200,50,800,600])
pause(2)
close all;













