function court_plot()

%% court
% 3pt line
X = linspace(-22, 22, 300);
Y = sqrt(22^2 + 8.95^2 - X.^2);

x_3ptline = [-22, -22, X, 22, 22] * 10;
y_3ptline = [0, 8.95, Y, 8.95, 0] * 9.85;

x_3ptline = [-220, x_3ptline, 220];
y_3ptline = [-50, y_3ptline, -50];

plot(x_3ptline, y_3ptline, '-k', 'LineWidth', 3); hold on;

% paint area
X = linspace(-6, 6, 300);
Y = sqrt(6^2 - X.^2);

x_circle = X * 10;
y_circle1 = (Y  + 15) * 9.85;
y_circle2 = (-Y  + 15) * 9.85;

plot(x_circle, y_circle1, '-k', 'LineWidth', 3);
plot(x_circle, y_circle2, '-k', 'LineWidth', 3);

x_paint1 = [-8, -8, 8, 8]*10;
y_paint1 = [-50, 15* 9.85, 15*9.85, -50];
plot(x_paint1, y_paint1, '-k', 'LineWidth', 3);

x_paint2 = [-6, -6, 6, 6]*10;
y_paint2 = [-50, 15* 9.85, 15*9.85, -50];
plot(x_paint2, y_paint2, '-k', 'LineWidth', 3);

% basket area
X = linspace(-0.75, 0.75, 300);
Y = sqrt(0.75^2 - X.^2);

x_basket = X * 10;
y_basket1 = (Y  + 0.75) * 9.85;
y_basket2 = (-Y  + 0.75) * 9.85;

plot(x_basket, y_basket1, '-k', 'LineWidth', 3);
plot(x_basket, y_basket2, '-k', 'LineWidth', 3);

% 3 sec area
X = linspace(-4, 4, 300);
Y = sqrt(4^2 - X.^2);

x_3sec = X * 10;
y_3sec = (Y + 0) * 9.85;
plot(x_3sec, y_3sec, '-k', 'LineWidth', 3); hold off;
set(gcf,'position',[200,50,800,600])
set(gca,'xtick',{}, 'ytick',{});




end