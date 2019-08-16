function out = objective_function_coefs_2(list_num_links, list_num_station_links, num_of_nodes, total_epochs, NUMBER_OF_SATELLITES)
    coefficients = null(1);
    for epoch = 1:total_epochs  
       num_flows = list_num_links(epoch)*2 - list_num_station_links(epoch); % afairw ta links satel to station
        coefficients = [ coefficients, ones(1, num_flows)*(10^(-3))];% flows
        coefficients = [ coefficients, zeros(1, NUMBER_OF_SATELLITES)]; % divergencies of satellites
        coefficients = [ coefficients, -ones(1, num_of_nodes - NUMBER_OF_SATELLITES)]; % divergencies of stations (should have negative coefficient in objective function)
%         coefficients = [ coefficients, zeros(1, num_of_nodes)]; % buffers <old version>
        coefficients = [ coefficients, ones(1, num_of_nodes)*(10^(-3))]; % buffers <new version> added delay minimization
    end
    out = coefficients;  
end