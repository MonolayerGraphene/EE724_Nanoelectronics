
% Define the range of x values
x = 4:0.05:16;

% Calculate corresponding y values for each function
y1 = 37.5*(0.2-(1./x));
y2 = 37.5./(x.^2);
y3 = 0.29*exp(-(x-15).^2 / (2 * 0.1^2)) / (0.1 * sqrt(2 * pi)) - 0.29*exp(-(x-5).^2 / (2 * 0.1^2)) / (0.1 * sqrt(2 * pi));

% Plot the functions
plot(x, y1, 'b-', x, y2, 'r-', x, y3, 'g-', 'LineWidth', 2);

% Add labels and title
xlabel('x');
ylabel('y');
title('Plot of V, E, and charge');

% Display the grid with increased detail
grid on;
grid minor; % Display a more detailed grid

% Add a legend
legend('Potential', 'Electric field', 'Charge');

% Display the plot
axis equal; % Optional: equal scaling for x and y axes

