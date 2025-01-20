% Define the range of x values
x = 1:0.1:16;

% Calculate corresponding y values for y = x
y = 37.5*(1./x-1/15);

% Plot the function
plot(x, y);

% Add labels and title
xlabel('r');
ylabel('E');
title('Energy variation with inner radius');

% Display the grid
grid on;
grid minor;

legend('Energy')
% Display the plot
axis equal; % Optional: equal scaling for x and y axes