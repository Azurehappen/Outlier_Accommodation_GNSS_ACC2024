function tidx = ephtidx(t_oc,t_sv,SV_health,message_duration)
% find the time index in eph data

dtr = limit_tgps(t_sv-t_oc); % Limit time (in seconds) to 1-week
tidx = find(dtr>=-message_duration&dtr<=message_duration);
% Satellite health check
i = SV_health(tidx)==0;
tidx = tidx(i);

end