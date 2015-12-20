function [ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param )

if isfield( param, 'cycleNUM' )
    cycleNUM = param.cycleNUM;
else
    cycleNUM = inf;
end

if cycleNUM > 0
    cycleNUM_S = cycleNUM;
    cycleNUM_E = cycleNUM;
else
    if isfield( param, 'cycleNUM_SAMPLE' )
        cycleNUM_S = param.cycleNUM_SAMPLE;
    else
        cycleNUM_S = inf;
    end
    if isfield( param, 'cycleNUM_EXEMPLAR' )
        cycleNUM_E = param.cycleNUM_EXEMPLAR;
    else
        cycleNUM_E = inf;
    end
end


end

