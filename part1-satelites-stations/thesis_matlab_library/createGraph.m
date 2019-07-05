%________ CREATE GRAPH FUNCTION ________
num_satellites = 2;
nodes = 1:3;
xij = {}; % Tha dinetai san paragraphos ths synarthshs
xij{1} = [ 1     1     2];
xij{2} = [-1     1     3];
xij{3} = [-1     2     3];
opt_results = [0.0004  100.0000  100.0000  0.0004 , 100.0000  100.0000  200.0000]; % edge weights
%%% ____________________________________

outgoing1 = null(1,1); outgoing2 = null(1,1);
incoming1 = null(1,1); incoming2 = null(1,1);
node_names = {}; % graph parameter
markers = {}; % graph parameter
unwanted = null(1,1);
for i = 1:length(xij)
    outgoing1 = [outgoing1, xij{i}(2)];
    incoming1 = [incoming1, xij{i}(3)];
    outgoing2 = [outgoing2, xij{i}(3)];
    incoming2 = [incoming2, xij{i}(2)];
    if xij{i}(1) == -1
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
% sanity-check:
disp('outgoing:') 
disp(outgoing)
disp('incoming:') 
disp(incoming)

G = digraph(outgoing , incoming, opt_results);
LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);
plot(G, 'EdgeLabel', G.Edges.Weight, 'EdgeColor', 'k', 'NodeLabel', node_names, 'Marker', markers, 'LineWidth', LWidths, 'MarkerSize',6,'NodeColor',[1 0 0])








