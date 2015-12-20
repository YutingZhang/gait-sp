function [param param_e idx] = zgParamEmutatorNext( param_e )

if param_e.empty_param > 0
    param_e.empty_param = -1;
    param = []; %empty
    idx = 1;
    return;
elseif param_e.empty_param <0
    param = {}; %EOF
    idx = [];
    return;
end

for k=1:(length(param_e.fields)-1)
    d = param_e.index(k)+1;
    while d > param_e.total(k)
        param_e.index(k+1) = param_e.index(k+1)+1;
        d = d - param_e.total(k);
    end
    param_e.index(k) = d - 1;
end

idx = param_e.index+1;

if ( param_e.index(end) >= param_e.total(end) )
    param = {}; %EOF
    return;
end

for k=1:length(param_e.fields)
    EXPR = ['param.' param_e.fields{k} ' = param_e.ranges.' param_e.fields{k} '(' int2str(param_e.index(k)+1) ');'];
    eval(EXPR);
end

param_e.index(1) = param_e.index(1) + 1;

end