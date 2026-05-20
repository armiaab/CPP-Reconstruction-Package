function OrganizeData(DOY,BGA,n_dailypoints)

% Initialize arrays to store the averaged quantity and the corresponding DOY
averaged_quantity = [];
averaged_DOY = [];
t1 = floor(DOY);

% Iterate over each day
for day = min(t1):max(t1)
    % Find the indices corresponding to the current day
    indices = find(DOY == day);
     
    % Calculate the number of data points per time point
    num_points_per_time_point = floor(numel(indices) / n_dailypoints);
    
    % Initialize arrays to store the averaged quantity and the corresponding DOY for the current day
    day_averaged_quantity = zeros(1, n_dailypoints);
    day_averaged_DOY = zeros(1, n_dailypoints);
    
    % Iterate over each time point
    n_dailypoints
    for t = 1:n_dailypoints
        % Calculate the starting and ending indices for the current time point
        start_index = (t - 1) * num_points_per_time_point + 1;
        end_index = min(start_index + num_points_per_time_point - 1, numel(indices));
        
        % Calculate the average quantity for the current time point
        day_averaged_quantity(t) = mean(BGA(indices(start_index:end_index)));
        
        % Calculate the corresponding DOY for the current time point
        day_averaged_DOY(t) = mean(DOY(indices(start_index:end_index)));
    end
    
    % Append the averaged quantity and the corresponding DOY for the current day to the overall arrays
    averaged_quantity = [averaged_quantity day_averaged_quantity];
    averaged_DOY = [averaged_DOY day_averaged_DOY];
end

% Plot the averaged quantity versus the corresponding DOY
size(averaged_DOY)
size(averaged_quantity)
nnz(isnan(averaged_DOY))
nnz(isnan(averaged_quantity))
plot(averaged_DOY, averaged_quantity, '.-');
xlabel('Day of Year');
ylabel('Averaged Quantity');
title('Averaged Quantity versus Day of Year');