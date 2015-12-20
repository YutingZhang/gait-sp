function feat_fun_param = zgInitParam_EigenStep( DATA, options )

authFoldN = 3;
eng_threshold = 0.95;
if exist( 'options', 'var' )
    if isfield( options, 'authFoldN' )
        authFoldN = options.authFoldN;
    end
end

FEAT = zgBatchExtractFeatures( DATA, @zgFeature_AllCycles );

person_total_both = size(FEAT.BOTH,2);
day_total    = size(FEAT.BOTH,3);

[ ~, authFoldRange ] = zgFoldPartition( person_total_both, authFoldN );

pca_model.ident_pca = cell(1,day_total);
pca_model.auth_pca  = cell(authFoldN,day_total);
for d=1:day_total
    pca_model.ident_pca{d} = BatchRecordPCA( FEAT.BOTH(:,:,d), eng_threshold );
    for k=1:authFoldN
        pca_model.auth_pca{k,d} = BatchRecordPCA( ...
            FEAT.BOTH(:,authFoldRange(k,1):authFoldRange(k,2),d), eng_threshold );
    end
end

feat_fun_param = { pca_model };

end

function p = BatchRecordPCA( RECORD, eng_threshold )

A  = cat( 1, RECORD{:} );
CN = size(A,2);
p  = cell(CN,1);
parfor c=1:CN
    X = cell2mat(A(:,c).');
    t = [];
    [t.Dic,t.Xm,t.Eng]  = zgPCASpace( X, eng_threshold );
    p{c} = t;
end

end

