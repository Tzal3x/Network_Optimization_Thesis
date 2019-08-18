function out = calculate_delay_v1(opt_divergencies_ca, NUMBER_OF_SATELLITES)
    % Calculate average total information received by stations per epoch.
    % Used to measure the delay in the network. We define delay as the 
    
    sum_stat_divs = null(1,1); % sum of station divergences PER EPOCH(total info received in given epoch)
    total_epochs = length(opt_divergencies_ca);
    num_nodes = length(opt_divergencies_ca{1});
    for i = 1:length(opt_divergencies_ca)
        stat_divs = opt_divergencies_ca{i}((NUMBER_OF_SATELLITES+1):num_nodes); % station divergences of given epoch
        sum_stat_divs = [sum_stat_divs, sum(stat_divs)];
    end
    res = sum(sum_stat_divs)/total_epochs;
    out = res;
end