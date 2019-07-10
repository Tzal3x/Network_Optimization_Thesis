% GraphOnMap(coords_ca, xij_ca, opt_results_ca, nodes)

% for i = 1:3
%    disp(i);
%    a = figure;
%    pause;
%    delete(a)
% end
coords_ca = {};% REMOVE THIS!!!
coords_ca{1} = [-44.2242, -244.6731,  -26.0628;
                -145.7537, -195.3200,  -55.7311;
                -138.7760, -134.8590,  -50.5399]; % coords of epoch 1        
total_epochs = length(coords_ca);
lons = null(1,1);
lats = null(1,1);
n = length(coords_ca{1}(:,1));
for i = 1:total_epochs
    for node = 1:n
        lons = [lons, coords_ca{i}(node,1)];
        lats = [lats, coords_ca{i}(node,2)];
    end
end
%% Create Map
    hold on
    m_proj('miller','lat',[-77 77]);   
    m_coast('patch',[.7 1 .7],'edgecolor','none'); 
    m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.2 .65 1]);

    cities={'Cairo','Washington','Buenos Aires'}; 
%     lons=[ 30+2/60  -77-2/60   -58-22/60];
%     lats=[ 31+21/60  38+53/60  -34-45/60];
    
    
    for k=1:3
      ln=lons(k);
      lt=lats(k);
      disp(ln);disp(lt);
%       m_line(ln,lt,'color','r','linewi',2); 
      m_patch(ln,lt)
%       m_text(ln(end),lt(end),sprintf('%s - %d km',cities{k},round(range)));
    end
    
    title('Great Circle Routes','fontsize',14,'fontweight','bold');
    
    set(gcf,'color','w');   % Need to do this otherwise 'print' turns the lakes black
    
    