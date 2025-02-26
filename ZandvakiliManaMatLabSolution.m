% Parameters
n_iterations = 50;  % Number of iterations
initial0 = 0.35;           % Initial value
bifurcation_values = [2.5, 3.2, 4.0];  % Values of bifurcation
bifurcation_range = linspace(2.5, 4, 1000); % Range of bifurcation values
% Storing and Set the initial value
l_fixed = zeros(1, n_iterations + 1);
l_periodic = zeros(1, n_iterations + 1);
l_chaotic = zeros(1, n_iterations + 1);
l_fixed(1) = initial0;
l_periodic(1) = initial0;
l_chaotic(1) = initial0;
% logistic map for the 2.5, 3.2 and 4 bifurcation values
for i = 1:length(bifurcation_values)
    r = bifurcation_values(i);
    l = initial0;
    for n = 2:n_iterations+1
        l = r * l * (1 - l);
        if r == 2.5
            l_fixed(n) = l;
        elseif r == 3.2
            l_periodic(n) = l;
        elseif r == 4.0
            l_chaotic(n) = l;
        end
    end
end
figure;
% Figure (a) 
subplot(2, 2, 1);
plot(0:n_iterations, l_fixed, 'ko-');  % Include the 0th iteration
title('r = 2.5 (Fixed Behavior)'); % Set title
xlabel('Iteration (n)'); % Set label for x
ylabel('l[n]'); % Set label for y
ylim([0 1]);
grid on;

% Figure (b) 
subplot(2, 2, 2);
plot(0:n_iterations, l_periodic, 'ko-');  % Include the 0th iteration
title('r = 3.2 (Periodic Behavior)');
xlabel('Iteration (n)');
ylabel('l[n]');
ylim([0 1]);
grid on;

% Figure (c) 
subplot(2, 2, 3);
plot(0:n_iterations, l_chaotic, 'ko-');  % Include the 0th iteration
title('r = 4 (Chaotic Behavior)');
xlabel('Iteration (n)');
ylabel('l[n]');
ylim([0 1]);
grid on;
% Parameters for bifurcation
scale = 700; % determines the level of rounding
max_values = 200; % maximum values to plot
interations = 500; % number of iterations of logistic equation
initial0 = 0.35; % initial value
% Loop through the "r" values and calculate logistic map
data = [];
for j = 1:length(bifurcation_range)
    r = bifurcation_range(j); % current "r" value
    l = zeros(interations, 1); % storing
    l(1) = initial0; % set initial condition
    
    % Iterate the logistic map
    for i = 2:interations
        l(i) = r * l(i-1) * (1 - l(i-1));
    end
    % Save unique stable values 
    unique_vals = unique(round(scale * l(end-max_values:end)));
    data = [data; bifurcation_range(j) * ones(length(unique_vals), 1), unique_vals];
end
% Create the subplot for Bifurcation Diagram in position (2, 2, 4)
subplot(2, 2, 4);
h = plot(data(:, 1), data(:, 2) / scale, 'k.');
set(h, 'markersize', 1);
title('Bifurcation Diagram');
xlabel('r');
ylabel('l[n]');
xlim([2.5 4]);
ylim([0 1]);
grid on;
