% Load the CSV file
data = readtable('joinville_transect_ctds_without_gamma_n.csv');

% Get unique events
uniqueEvents = unique(data.Event);

% Initialize an empty table to store results
resultTable = data;

% Loop through each unique event, which corresponds to a profile 
for i = 1:length(uniqueEvents)
    % Get the current event
    currentEvent = uniqueEvents{i};
    % Print the current event
    disp(['Processing event: ', currentEvent]);
    
    % Find rows with the current event
    eventRows = strcmp(data.Event, currentEvent);

    % Apply eos80_legacy_gamma_n to the current profile
    for j = find(eventRows)'
        % Extract the necessary profile data
        profile_in_situ_temp = data.Temp__C_(j);        % in-situ temperature (ITS-90) [ deg C ]
        profile_pressure = data.Press_dbar_(j);         % sea pressure [ dbar ]   
        profile_practical_salinity = data.Sal(j);       % Practical Salinity [ unitless ]
        profile_lat = data.Latitude(j);
        profile_lon = data.Longitude(j);

        % Apply the function (modify this according to your actual function)
        % order of parameters (SP,t,p,long,lat)
        profile_neutral_density = eos80_legacy_gamma_n(profile_practical_salinity, profile_in_situ_temp, profile_pressure, profile_lon, profile_lat);       
        %[profile_neutral_density, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n(profile_practical_salinity, profile_in_situ_temp, profile_pressure, profile_lon, profile_lat);       
        % Store the result in a new column (or modify the existing column)
        resultTable.neutral_density{j} = profile_neutral_density;
    end
end

% Save the result back to a CSV file
writetable(resultTable, 'joinville_transect_CTDs_neutral_density.csv');
disp('done!');