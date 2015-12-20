function P = zgStr2Param( s )

SPLITTER_IDX = find(s == ';');
AUG_SPLITTER_IDX = [ (diff( [0 SPLITTER_IDX] ))-1 ; ones(1,length(SPLITTER_IDX)) ];
STR_CELL = mat2cell( s(1:SPLITTER_IDX(end)), 1, reshape( AUG_SPLITTER_IDX, numel(AUG_SPLITTER_IDX), 1) ) ;
B = reshape( STR_CELL, 2, length(SPLITTER_IDX) );
B = B(1,:);
P = [];

if length(B)==1 && strcmp( B{1}, 'default' )
    P = [];
    return;
end

for k=1:length(B)
    c = B{k};
    eqidx = find(c=='-',1);
    if isempty(eqidx)
        error('no value is specified');
    end
    fname  = c(1:eqidx-1);
    fvalue = str2num( c(eqidx+1:end) );
    eval( [ 'P.' fname '= fvalue;' ] );
end


end
