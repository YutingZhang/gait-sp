function ps = zgParam2Str( param )
% e.g. param.x = 1:10, param.y = [2 4]

if isempty( param )
    ps = 'default;';
    return;
end

fields = fieldnames( param );
ps = '';
for k=1:length(fields)
    eval(['v = param.' fields{k} ';']);
    ps = [ ps fields{k} '-' num2str(v) ';' ];
end

end