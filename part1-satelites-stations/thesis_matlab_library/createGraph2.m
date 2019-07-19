function createGraph2(num_satellites, xij_ca, opt_results_ca, distances_ca, opt_buffers_ca, nodes, epoch_name)
% WARNING, LEFT UNDER CONSTRUCTION 17-7-2019. DO NOT USE IT!
% Version 2, added info about distances and buffers
outgoing1 = null(1,1); outgoing2 = null(1,1);
incoming1 = null(1,1); incoming2 = null(1,1);
node_names = {}; % graph parameter
markers = {}; % graph parameter
unwanted = null(1,1);
dists = null(1,1);
for i = 1:length(xij_ca)
    p1 = xij_ca{i}(2); % parent 1
    p2 = xij_ca{i}(3); % parent 2
    outgoing1 = [outgoing1, p1];
    incoming1 = [incoming1, p2];
    outgoing2 = [outgoing2, p2];
    incoming2 = [incoming2, p1];
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

for i = 1:length(outgoing)
    if i > length([outgoing1, outgoing2])
         dists = [dists, ' '];
    else
        dists = [dists, string(distances_ca(outgoing(i), incoming(i)))];
    end
end
dists(1:length([outgoing1, outgoing2])) = " \ dist: " + dists(1:length([outgoing1, outgoing2]));
% sanity-check:
disp('opt_results_ca:')
disp(opt_results_ca)
disp('outgoing:') 
disp(outgoing)
disp('incoming:') 
disp(incoming)
disp('dists:')
disp(dists)

% Create node characteristics: - - - - - - - - - - - 
% buf_history = {}; % buffer value of previous and current epoch for each node
% for i = 1:length(nodes)
%     if epoch_name ~= 0
%       buf_history{i} = ['s('+string(epoch_name)+'-1):'+opt_buffers_ca{i} + newline + ];
%     else 
%         
%     end
% end

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

G = digraph(outgoing , incoming, opt_results_ca);
LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);
disp('G.Edges.Weight: (PROBLEM ! WRONG ORDER)')% debug
disp(G.Edges.Weight) %debug
edgelabels = cellstr(G.Edges.Weight + dists');
disp('edgelabels:') % debug
disp(edgelabels) % debug

f = figure; %#ok<NASGU>
plot(G, 'EdgeLabel', edgelabels, 'EdgeColor', 'k', 'NodeLabel', node_names,...
    'Marker', markers, 'LineWidth', LWidths, 'MarkerSize',6,'NodeColor',[1 0 0],...
    'ArrowSize',17)

if epoch_name ~= 0
    title('DTNUM Graph (epoch '+string(epoch_name)+')')
else
    title('DTNUM Graph')
end

end