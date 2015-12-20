function FEAT = zgBatchExtractFeatures( DATA, feature_fun, ff_options, feature_fun_param ) 

disp('Extracting features ...');
N = numel(DATA.BOTH) + numel(DATA.SINGLE);
pb = zgProgressBarInit(N);

if ~exist( 'feature_fun_param', 'var' )
    feature_fun_param = {};
elseif ~iscell( feature_fun_param )
    feature_fun_param = { feature_fun_param };
end

recordNum_fun = @zgSingleRecordNum;
if exist( 'ff_options', 'var' )
    if isfield( ff_options, 'useCyclicAllRecord' )
        if ff_options.useCyclicAllRecord
            recordNum_fun = @zgCyclicRecordNum;
        end
    end
end

record_total = size(DATA.BOTH,1);
for k3=1:size(DATA.BOTH,3)
    for k2=1:size(DATA.BOTH,2)
        for k1=1:record_total
            K1 = recordNum_fun(k1,record_total);
            FEAT.BOTH{k1,k2,k3} = feature_fun( cell2mat(DATA.BOTH(K1,k2,k3)), feature_fun_param{:} );
        end
        pb = zgProgressBarNext( pb,size(DATA.BOTH,1) );
    end
end

record_total = size(DATA.SINGLE,1);
for k2=1:size(DATA.SINGLE,2)
    for k1=1:record_total
        K1 = recordNum_fun(k1,record_total); 
        FEAT.SINGLE{k1,k2} = feature_fun( cell2mat(DATA.SINGLE(K1,k2)), feature_fun_param{:} );
    end
    pb = zgProgressBarNext( pb,size(DATA.SINGLE,1) );
end


end

function r = zgSingleRecordNum(r, total)
end

function c = zgCyclicRecordNum(r, total)
c=[r:total 1:(r-1)];
end
