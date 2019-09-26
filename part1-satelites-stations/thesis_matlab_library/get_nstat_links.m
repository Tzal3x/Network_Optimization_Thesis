function out = get_nstat_links(xij, NUMBER_OF_SATELLITES)
% Returns the number of links that have a station as a parent
    num_station_links = 0; % number of links that have a station as a parent
    for i = 1:length(xij)
        if xij{i}(2) > NUMBER_OF_SATELLITES || xij{i}(3) > NUMBER_OF_SATELLITES
            num_station_links = num_station_links + 1;
        end
    end
    out = num_station_links;
end