%% --------- HELPING SCRIPT FOR EXPERIMENTATION AND DEBUGING. IGNORE IT -----------------------

load coastlines
axesm mollweid
framem('FEdgeColor','blue','FLineWidth',0.5)

a = plotm(coastlat,coastlon,'LineWidth',1,'Color','blue');

citylats = [30 -23 -32]; citylongs = [32 -43 116];
b = plotm(citylats(1),citylongs(1),'r.','MarkerSize',15);
c = plotm(citylats(2),citylongs(2),'r.','MarkerSize',15);
% d = plotm(citylats(3),citylongs(3),'r.','MarkerSize',15);

% Connect lines
% [gclat,gclong] = track2('gc',citylats(1),citylongs(1),...
%                              citylats(2),citylongs(2)); %gc means great circle track
%                          
% [rhlat,rhlong] = track2('rh',citylats(1),citylongs(1),...
%                              citylats(3),citylongs(3));% rh: rhumb line track
%                          
% plotm(gclat,gclong,'m-'); plotm(rhlat,rhlong,'m-')


% con = plot([citylats(1) citylongs(1)],[citylats(3)+i citylongs(3)],'r');
% % connecting lines this way will not work because it uses the cartesian
% coordinate system

for i=1:100
    d = plotm(citylats(3)+i,citylongs(3),'r.','MarkerSize',20);
    [rhlat,rhlong] = track2('rh',citylats(2),citylongs(2),citylats(3)+i,citylongs(3));% rh: rhumb line track
    conn = plotm(rhlat,rhlong,'m-');
    pause(1)
    delete(d)
    delete(conn)
end

