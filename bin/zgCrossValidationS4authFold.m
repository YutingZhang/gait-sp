function [R W paramPersonNum] = zgCrossValidationS4authFold( FEAT , score_fun , param, options )
% R = zgCrossValidationS( FEAT , score_fun , param )
if ischar(FEAT)
    load( FEAT );
end

paramPersonNum.both  = size(FEAT.BOTH,2);
paramPersonNum.single = size(FEAT.SINGLE,2);

%% Override score_fun
% score_fun = @zgScoreArr_SPVoting;

%% Init num
person_total_both   = size(FEAT.BOTH,2);
record_total = size(FEAT.BOTH,1);
person_total_single = size(FEAT.SINGLE,2);
person_total = person_total_both+person_total_single;

%% Parsing options
authFoldN = 3;
norm_method = 'none';
if exist( 'options', 'var' )
    if isfield( options, 'authFoldN' )
        authFoldN = options.authFoldN;
    end
    if isfield( options, 'PreNormFunName' )
        norm_method = options.PreNormFunName;
    end
end
[ authFoldEleN, authFoldRange ] = zgFoldPartition( person_total_both, authFoldN );


%% Process

% channel, sample record number, exemplar record number,
% sample persons, exemplar persons, days
R.IDENT = zeros( [ 5, record_total , record_total , ...
    person_total_both , person_total, 2 ] );
ts = size(R.IDENT);
channelN = ts(1);
W.IDENT  = zeros(ts);

R.AUTH   = cell(authFoldN,1); W.AUTH   = cell(authFoldN,1);
R.AUTH_NORM = cell(authFoldN,1); W.AUTH_NORM = cell(authFoldN,1);
cohortIDXArr = false(person_total,authFoldN);
authFoldSamN = zeros(authFoldN,1);
authFoldExeN = zeros(authFoldN,1);
for k = 1:authFoldN
    authFoldSamN(k) = person_total-authFoldEleN(k);
    authFoldExeN(k) = person_total_both-authFoldEleN(k);
    R.AUTH{k} = zeros( [ 5, record_total , record_total , ...
        authFoldSamN(k) , person_total_both-authFoldEleN(k), 2 ] );
    fts = size(R.AUTH{k});
    W.AUTH{k} = zeros(fts);
    R.AUTH_NORM{k} = zeros(fts);
%     W.AUTH_NORM{k} = zeros(fts);
    cohortIDXNum = authFoldRange(k,1):authFoldRange(k,2);
    cohortIDXArr( cohortIDXNum, k ) = true;
end
samIDXArr      = ~cohortIDXArr;
exeIDXArr = samIDXArr;
exeIDXArr( (person_total_both+1):end, : ) = false;
exeIDXNum = sum(exeIDXArr,1);


R1 = zeros( [ 5, 1 , 1 , ...
    person_total , person_total ] );
W1 = zeros( size(R1) );
for d1 = 1:2
    d2 = 3-d1;
    param.cvsDayID = d1;    % exemplarDay
    for er=1:record_total
        E_BOTH   = FEAT.BOTH(er,:,d1);
        E_SINGLE = FEAT.SINGLE(er,:);
        for sr=1:record_total
            S_BOTH = FEAT.BOTH(sr,:,d2);
            S_SINGLE = FEAT.SINGLE(sr,:);
            %% identification
            param.cvsRecType = 'ident';
            parfor ek = 1:person_total_both
                R0 = zeros([channelN,1,1,person_total]);
                W0 = zeros( size(R0) );
                exemplar = E_BOTH{1,ek};
                samples  = S_BOTH;
                for sk=1:person_total_both
                    sample = samples{1,sk};
                    [ R0( :, 1, 1, sk ), W0( :, 1, 1, sk ) ] = feval( score_fun, sample, exemplar, param );
                end
                R1(:,1,1,:, ek) = R0;
                W1(:,1,1,:, ek) = W0;
            end
            parfor ek = 1:person_total_single
                R0 = zeros([channelN,1,1,person_total]);
                W0 = zeros( size(R0) );
                exemplar = E_SINGLE{1,ek};
                samples  = S_BOTH;
                for sk=1:person_total_both
                    sample = samples{1,sk};
                    [ R0( :, 1, 1, sk ), W0( :, 1, 1, sk ) ] = feval( score_fun, sample, exemplar, param );
                end
                R1(:, 1,1,:, ek+person_total_both) = R0;
                W1(:, 1,1,:, ek+person_total_both) = W0;
            end
            R.IDENT( :, sr, er, :, :, d1 ) = R1(:,1,1,1:person_total_both,:);
            W.IDENT( :, sr, er, :, :, d1 ) = W1(:,1,1,1:person_total_both,:);
            
            %% authentication
            param.cvsRecType = 'auth';
            for k = 1:authFoldN
                param.cvsFoldID  = k;
                parfor ek = 1:person_total_both
                    R0 = zeros([channelN,1,1,person_total]);
                    W0 = zeros( size(R0) );
                    exemplar = E_BOTH{1,ek};
                    samples = S_BOTH;
                    for sk=1:person_total_both
                        sample = samples{1,sk};
                        [ R0( :, 1, 1, sk ), W0( :, 1, 1, sk ) ] = feval( score_fun, sample, exemplar, param );
                    end
                    samples  = S_SINGLE;
                    for sk=1:person_total_single
                        sample = samples{1,sk};
                        [ R0( :, 1, 1, sk+person_total_both ), W0( :, 1, 1, sk+person_total_both ) ] = feval( score_fun, sample, exemplar, param );
                    end
                    R1(:,1,1,:, ek) = R0;
                    W1(:,1,1,:, ek) = W0;
                end
                R.AUTH{k}( :, sr, er, :, :, d1 ) = R1(:,1,1,samIDXArr(:,k),exeIDXArr(:,k));
                W.AUTH{k}( :, sr, er, :, :, d1 ) = W1(:,1,1,samIDXArr(:,k),exeIDXArr(:,k));
%                 W.AUTH_NORM{k}( :, sr, er, :, :, d1 ) = W.AUTH{k}( :, sr, er, :, :, d1 );
                [NormCallBack R1N] = zgScoreNormInit( R1(:,1,1,:,1:person_total_both), norm_method );
                R.AUTH_NORM{k}( :, sr, er, :, :, d1 )  = feval( NormCallBack, ...
                    exp( log( R1N(:,1,1,samIDXArr(:,k),exeIDXArr(:,k)) ) - ...
                    log(repmat( mean(R1N(:,1,1,samIDXArr(:,k),exeIDXArr(:,k)),5), ...
                    [ 1 1 1 1 exeIDXNum(k) ] )) ) ...
                );
            end

        end
        fprintf(1,'.');
    end
end

for k = 1:authFoldN
    W.AUTH_NORM{k} = W.AUTH{k};
end

fprintf(1,'\n');

end

