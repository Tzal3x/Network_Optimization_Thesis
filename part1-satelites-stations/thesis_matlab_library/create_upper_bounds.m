function out = create_upper_bounds(list_num_links, list_num_station_links, num_of_nodes, total_epochs)
% Creates upper_bound vector for problem of "part_2"
    bounds = null(1);
    flow_capacity = 20;
    buffer_capacity = 10;
    for epoch = 1:total_epochs      
        num_flows = list_num_links(epoch)*2 - list_num_station_links(epoch);
        bounds = [ bounds, zeros(1, num_flows) + flow_capacity];% flows
        bounds = [ bounds, zeros(1, num_of_nodes) + inf]; % divergencies
        bounds = [ bounds, zeros(1, num_of_nodes) + buffer_capacity]; % buffers
    end
    out = bounds;   
end