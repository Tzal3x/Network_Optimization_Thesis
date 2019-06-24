function out = create_upper_bounds(num_flows, num_of_nodes, total_epochs)
% Creates upper_bound vector for problem of "part_2"
    bounds = null(1);
    flow_capacity = 20;
    buffer_capacity = 10;
    for epoch = 1:total_epochs      
        bounds = [ bounds, zeros(1, num_flows) + flow_capacity];% flows
        bounds = [ bounds, zeros(1, num_of_nodes) + inf]; % divergencies
        bounds = [ bounds, zeros(1, num_of_nodes) + buffer_capacity]; % buffers
    end
    out = bounds;   
end