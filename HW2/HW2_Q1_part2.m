% Define the range of x values
x1 = 4:0.1:10;
x2 = 10:0.1:16;

% Define the functions for the first interval

y11 = 45*(0.2-(1./x1));
y12 = 15*(0.1-(1./x2))+4.5;         
        % y = x^2 for 1<=x<=5

% Define the functions for the second interval
y21 = 45./(x1.^2);        % y = x^2 for 5<=x<=10
y22 = 15./(x2.^2);        % y = x^3 for 5<=x<=10

% Define the functions for both intervals
y31 = exp(-(x1-5).^2 / (2 * 0.1^2)) / (0.1 * sqrt(2 * pi))-0.8*exp(-(x1-10).^2 / (2 * 0.1^2)) / (0.1 * sqrt(2 * pi));           % y = x for 1<=x<=5
y32 = -0.8*exp(-(x2-15).^2 / (2 * 0.1^2)) / (0.1 * sqrt(2 * pi));      % y = log(x) for 1<=x<=5

% Combine the pieces
x = [x1, x2];
ya = [y11, y12];
yb = [y21, y22];
yc = [y31, y32];

% Plot the piecewise functions
hold on
plot(x, ya, 'LineWidth', 2);
plot(x, yb, 'LineWidth', 2);
plot(x, yc, 'LineWidth', 2);

% Add labels and title
xlabel('x');
ylabel('y');
title('Plot of V, E, and charge');

% Display the grid
grid on;
grid minor;
legend('Potential', 'Electric field', 'Charge');


% Display the plot

