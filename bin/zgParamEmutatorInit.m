function param_e = zgParamEmutatorInit( param_ranges )
% e.g. param.x = 1:10, param.y = [2 4]

if isempty(param_ranges)
    param_e.empty_param = 1;
    param_e.ranges = [];
    param_e.fields = {};
    param_e.index  = 0;
    param_e.total  = 1;
    return;
end
param_e.empty_param = 0;

param_e.ranges = param_ranges;
param_e.fields = fieldnames( param_ranges );
param_e.index  = zeros( 1, length(param_e.fields) );
param_e.total  = zeros( 1, length(param_e.fields) );
for k=1:length(param_e.fields)
    eval( ['L = length( param_e.ranges.' param_e.fields{k} ' );'] );
    if L<1
        error(['The range of ' param_e.fields{k} 'is empty']);
    end
    param_e.total(k) = L;
end

end