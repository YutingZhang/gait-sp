function D = zgFeature_EigenStep( RECORD, pca_model )

NormalizedLen = 100;
maxCycleNum   = 10;

CN  = size(RECORD.data,3);
CYC = zgFeature_AllCycles( RECORD, NormalizedLen );

cycleNumList = [1:maxCycleNum inf];

D = cell(CN,1);

for c=1:CN
    cyc = CYC{c};
    D{c}.cycs = cyc;
    D{c}.Y = cell( length(cycleNumList), 1 );
    for k=1:length(cycleNumList)
        n = min( size(cyc,2), cycleNumList(k) );
        X = mean(cyc(:,1:n), 2);
        D{c}.Y{k}.ident = zgCellFun( @transPCA, pca_model.ident_pca, {X,c} );
        D{c}.Y{k}.auth  = zgCellFun( @transPCA, pca_model.auth_pca, {X,c} );
    end
end


end


function Y = transPCA( p, X, channelID )

Y.y = zgPCATrans( X, p{channelID}.Dic, p{channelID}.Xm );
Y.eng = p{channelID}.Eng;

end

function R = zgCellFun( fun_handle , C, param_cell )

R = cell( size(C) );
n = numel(C);
for k=1:n
    R{k} = feval( fun_handle, C{k}, param_cell{:} );
end

end
