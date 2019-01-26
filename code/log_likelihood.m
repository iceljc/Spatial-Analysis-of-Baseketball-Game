function L = log_likelihood(zv)

z0 = 0.0;
% N = 50;
% Xv = zeros(1, (N-1)*(N-1));
% Xv = cell2mat(varargin);

%%

% f = xlsread('nba_savant201939.csv');
% 
% % label = f(:, 8);
% x = f(:, 13);
% y = f(:, 14);
% shot_distance = f(:, 11);
% 
% 
% index_long_dist = find(shot_distance >= 35);
% % label(index_long_dist) = [];
% x(index_long_dist) = [];
% y(index_long_dist) = [];
% 
% % index_shotmade = find(label == 1);
% % index_shotmiss = find(label == 0);
% 
% % x_made = x(index_shotmade);
% % y_made = y(index_shotmade);
% % 
% % x_miss = x(index_shotmiss);
% % y_miss = y(index_shotmiss);
%     
% N = 50;
% xd = linspace(-250, 250, N);
% yd = linspace(-50, 350, N);
% 
% xq = x;
% yq = y;
% 
% count = zeros(N-1, N-1);
% 
% for j = 1 : N-1
%     for i = 1 : N-1
%         xv = [xd(i), xd(i+1), xd(i+1), xd(i), xd(i)];
%         yv = [yd(j), yd(j), yd(j+1), yd(j+1), yd(j)];
%         in = inpolygon(xq, yq, xv, yv);
%         count(i,j) = sum(in(:));
% 
%     end
% end
% dA = (xv(2)-xv(1))*(yv(3)-yv(1));
% count = rot90(count);
% 
% 
% Xv = reshape(count', 1, (N-1)*(N-1));

% z0 = mean(Xv); 
dA = 1;

load main

%%

log_like = zeros(1, size(Xv, 2));


for i = 1 : size(Xv, 2)
    
    log_like(i) = -dA * exp(zv(i) + z0) + Xv(i) * (zv(i) + z0) + Xv(i) * log(dA) - log(factorial(Xv(i)));
    
end


L = sum(log_like);

end