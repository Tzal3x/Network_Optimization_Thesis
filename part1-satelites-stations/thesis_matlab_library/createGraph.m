function createGraph(num_satellites, xij_ca, opt_results_ca, nodes, epoch_name)
% Version 1
outgoing1 = null(1,1); outgoing2 = null(1,1);
incoming1 = null(1,1); incoming2 = null(1,1);
node_names = {}; % graph parameter
markers = {}; % graph parameter
unwanted = null(1,1);
for i = 1:length(xij_ca)
    outgoing1 = [outgoing1, xij_ca{i}(2)];
    incoming1 = [incoming1, xij_ca{i}(3)];
    outgoing2 = [outgoing2, xij_ca{i}(3)];
    incoming2 = [incoming2, xij_ca{i}(2)];
    if xij_ca{i}(1) == -1
       unwanted = [unwanted, i]; 
    end
end
% Remove station to satellite edges
outgoing2(unwanted) = [];
incoming2(unwanted) = [];
divs1 = [(length(nodes)+1):(length(nodes)+num_satellites),    (num_satellites+1):length(nodes)]; %source
divs2 = [1:num_satellites,   (length(nodes) + 1 + num_satellites):(length(nodes) * 2)]; %sink
outgoing = [outgoing1, outgoing2, divs1];
incoming = [incoming1, incoming2, divs2];

% Create node characteristics:
for i = 1:(length(nodes)*2) % *2 to add divergencies
    if i <= num_satellites
        node_names{i} = 'satellite_'+string(i);
        markers{i} = 'p';
    elseif i > num_satellites && i <= length(nodes)
        node_names{i} = 'station_'+string(i);
        markers{i} = 'o';  
    else
        node_names{i} = 'external_node_'+string(i);
        markers{i} = 'v';
    end
end

node_names = cellstr(node_names);
% % sanity-check:
% disp('outgoing:') 
% disp(outgoing)
% disp('incoming:') 
% disp(incoming)

G = digraph(outgoing , incoming, opt_results_ca);
LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);
f = figure;
plot(G, 'EdgeLabel', G.Edges.Weight, 'EdgeColor', 'k', 'NodeLabel', node_names,...
    'Marker', markers, 'LineWidth', LWidths, 'MarkerSize',6,'NodeColor',[1 0 0],...
    'ArrowSize',17)

if epoch_name ~= 0
    title('DTNUM Graph (epoch '+string(epoch_name)+')')
else
    title('DTNUM Graph')
end

end