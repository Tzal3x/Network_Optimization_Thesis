classdef station
    %Class describing a station object. Every station is a sink in graph
    %terms.
    properties
        longitude = 0;
        latitude = 0;
        neighbours = [] %binary list (of size N+M) containing the other sources (or sinks) that have a connection with the station
    end
    methods
        function station_obj = station(arg_long, arg_lat)
            station_obj.longitude = arg_long;
            station_obj.latitude = arg_lat;
        end
    end
end