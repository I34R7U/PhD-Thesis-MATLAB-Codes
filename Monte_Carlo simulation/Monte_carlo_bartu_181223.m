p = 9;
num_trials = 5e6;
positions = zeros(num_trials,p);
overlap = zeros(num_trials,1);
no_overlap = overlap;
P_overlap = 0.15;
for i=1:num_trials
    positions(i,1)=1;
    overlap_index = 1;
    for j=2:p
        random_number = rand;
        if random_number<=P_overlap
            positions(i,overlap_index) =  positions(i,overlap_index)+1;
        else
            positions(i,j) =  positions(i,j)+1;
            overlap_index =j;
        end
    end
end
for i=1:size(positions,1)
    for j=1:size(positions,2)
        if positions(i,j)>1
            overlap(i) = overlap(i)+1;
        end
    end
end        
% Plot histogram
histogram(overlap, 'Normalization', 'pdf'); % 'pdf' for probability density function
% Hold on to the current axes
hold on
% Fit a normal distribution to the data
%mu = mean(overlap);
%sigma = std(overlap);
%x = linspace(min(overlap), max(overlap), 100);
%y = normpdf(x, mu, sigma);
% Plot the normal fit
%plot(x, y, 'r', 'LineWidth', 2);
% Find the mode (most probable value)
%[~, mode_index] = max(y);
%mode_value = x(mode_index);
% Display the most probable value on the graph
%text(mode_value, max(y)*1.1, sprintf('Most probable: %.2f', mode_value), ...
 %   'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
  %  'BackgroundColor', 'w', 'EdgeColor', 'k');
% Draw a vertical line at the mode
%line([mode_value, mode_value], [0, max(y)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
% Set y-axis limits to [0, 1]
ylim([0, 1]);
% Add legend
%legend('Histogram', 'Gaussian Fit');
% Release the hold on the axes
hold off
% Add labels and title
xlabel('Number of sparks');
ylabel('Probability Density');
title('21 Pulses');

