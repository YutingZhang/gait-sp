function [R W paramPersonNum] = zgCrossValidationS4SPsrc( FEAT , score_fun , param, options )
% R = zgCrossValidationS( FEAT , score_fun , param )
if ischar(FEAT)
    load( FEAT );
end

paramPersonNum.both   = size(FEAT.BOTH,2);
paramPersonNum.single = size(FEAT.SINGLE,2);

channelN = length(FEAT.BOTH{1});

[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

%% Override score_fun
% score_fun = @zgScoreArr_SPVoting;

%% Init num
person_total_both = size(FEAT.BOTH,2);
record_total = size(FEAT.BOTH,1);
person_total_single = size(FEAT.SINGLE,2);
person_total = person_total_both+person_total_single;

%% Parsing options
authFoldN     = 3;
randSeed      = 0;
wordNumBase   = 10;
chosenClusterNum  = 3;
chosenExemplarNum = 30;
clusterCodingTol   = 0.1;
exemplarCodingTol  = 0.01;
softmaxSigma2      = 1;
hasIdent = 1;
hasAuth  = 1;
srcType  = 1;
listER   = 1:record_total;
detailedOutputPath = '';
if exist( 'options', 'var' )
    if isfield( options, 'authFoldN' )
        authFoldN = options.authFoldN;
    end
    if isfield( options, 'randSeed' )
        randSeed = options.randSeed;
    end
    if isfield( options, 'srcType' )
        srcType = options.srcType;
    end
    if isfield( options, 'detailedOutputPath' )
        detailedOutputPath = options.detailedOutputPath;
    end
    if isfield( options, 'listER' )
        listER  = options.listER;
    end
end
if ~isempty( detailedOutputPath )
    if ~exist( detailedOutputPath, 'dir' )
        mkdir( detailedOutputPath );
    end
end
if exist( 'param', 'var' )
    if isfield( param, 'wordNumBase' )
        wordNumBase = param.wordNumBase;
    end
    
    if isfield( param, 'chosenClusterNum' )
        chosenClusterNum = param.chosenClusterNum;
    end
    if isfield( param, 'chosenExemplarNum' )
        chosenExemplarNum = param.chosenExemplarNum;
    end

    if isfield( param, 'clusterCodingTol' )
        clusterCodingTol  = param.clusterCodingTol;
    end
    if isfield( param, 'exemplarCodingTol' )
        exemplarCodingTol = param.exemplarCodingTol;
    end
    if isfield( param, 'softmaxSigma' )
        softmaxSigma2      = (param.softmaxSigma).^2;
    end
    
    if isfield( param, 'hasIdent' )
        hasIdent = param.hasIdent;
    end
    if isfield( param, 'hasAuth' )
        hasAuth  = param.hasAuth;
    end
    
end
[ authFoldEleN, authFoldRange ] = zgFoldPartition( person_total_both, authFoldN );


%% Initialize

% ident
% channel, sample record number, exemplar record number,
% sample persons, exemplar persons, days
R.IDENT = zeros( [ channelN, record_total , record_total , ...
    person_total_both , person_total, 2 ] );
W.IDENT  = ones( size(R.IDENT) );

    
% auth
cohortIdxB   = false(person_total,authFoldN);
cohortIdx    = cell(authFoldN,1);
samIdx       = cell(authFoldN,1);
exeIdx       = cell(authFoldN,1);

authFoldSamN = person_total      - authFoldEleN;
authFoldExeN = person_total_both - authFoldEleN;

personFoldIdx   = zeros(1,person_total);
personIdxInFold = zeros( authFoldN, person_total );

for k = 1:authFoldN
    personFoldIdx(authFoldRange(k,1):authFoldRange(k,2)) = k;

    cohortIdx{k}   = authFoldRange(k,1):authFoldRange(k,2);
    cohortIdxB( cohortIdx{k}, k ) = true;
    samIdx{k}      = find( ~cohortIdxB( :, k ) ).';
    exeIdx{k}      = find( ~cohortIdxB( 1:person_total_both, k ) ).';
    personIdxInFold( k, samIdx{k} ) = 1:length(samIdx{k});
end

if hasAuth
    R.AUTH  = cell(authFoldN,1); W.AUTH = cell(authFoldN,1);
    for k = 1:authFoldN
        fts = [ channelN, record_total , record_total , ...
            authFoldEleN(k)+1, authFoldSamN(k) , authFoldExeN(k), 2 ];

        R.AUTH{k} = zeros( fts );
        W.AUTH{k} = ones( fts );
    end
end
    


%% Process


switch srcType
    case 1
        tolParameter = chosenExemplarNum;
    case 2
        tolParameter = exemplarCodingTol;
    case 3
        tolParameter = chosenExemplarNum;
    otherwise
        error( 'Unrecognized srcType' );
end

% sp portion
originalSPNum = 0;
selectedSPNum = 0;
for d1 = 1:2
    for r=1:record_total
        for c = 1:channelN
            for k=1:person_total_both
                originalSPNum = originalSPNum + length( FEAT.BOTH{r,k,d1}{c}.cycleID );
                selectedSPNum = selectedSPNum + sum( FEAT.BOTH{r,k,d1}{c}.cycleID <= cycleNUM_E );
            end
            for k=1:person_total_single
                originalSPNum = originalSPNum + length( FEAT.SINGLE{r,k}{c}.cycleID );
                selectedSPNum = selectedSPNum + sum( FEAT.SINGLE{r,k}{c}.cycleID <= cycleNUM_E );
            end
        end
    end
end
portionSP = selectedSPNum/originalSPNum;

spNum_s = zeros( channelN, record_total, person_total, 2 );
for c = 1:channelN
    % stack the data
    Ae = cell(record_total,2); % dictionary
    Be = cell(record_total,2); % label
    As = cell(record_total,2); % dictionary
    Bs = cell(record_total,2); % label
    Le = cell(record_total,2); % location for detailed analysis only
    
    foldExeIdxB   = cell(authFoldN,record_total,2);
    foldEleIdxB   = cell(authFoldN,record_total,2);
    foldEleIdx    = cell(authFoldN,record_total,2); % fold
    foldEleNum    = zeros(authFoldN,record_total,2);
    for d1 = 1:2
        for r=1:record_total
            Ae{r,d1} = cell(person_total,1); % dictionary
            Be{r,d1} = cell(person_total,1); % label
            As{r,d1} = cell(person_total,1); % dictionary
            Bs{r,d1} = cell(person_total,1); % label
            
            Le{r,d1} = cell(person_total,1); % location for detailed analysis only
            for k=1:person_total_both
                Ae{r,d1}{k} = FEAT.BOTH{r,k,d1}{c}.descriptor( FEAT.BOTH{r,k,d1}{c}.cycleID<=cycleNUM_E, : );
                Be{r,d1}{k} = repmat( k, size(Ae{r,d1}{k},1), 1 );
                As{r,d1}{k} = FEAT.BOTH{r,k,d1}{c}.descriptor( FEAT.BOTH{r,k,d1}{c}.cycleID<=cycleNUM_S, : );
                Bs{r,d1}{k} = repmat( k, size(As{r,d1}{k},1), 1 );
                
                Le{r,d1}{k} = FEAT.BOTH{r,k,d1}{c}.loc( FEAT.BOTH{r,k,d1}{c}.cycleID<=cycleNUM_E );
            end
            for k=1:person_total_single
                k2 = k + person_total_both;
                Ae{r,d1}{k2} = FEAT.SINGLE{r,k}{c}.descriptor( FEAT.SINGLE{r,k}{c}.cycleID<=cycleNUM_E, : );
                Be{r,d1}{k2} = repmat( k2, size(Ae{r,d1}{k2},1), 1 );
                As{r,d1}{k2} = FEAT.SINGLE{r,k}{c}.descriptor( FEAT.SINGLE{r,k}{c}.cycleID<=cycleNUM_S, : );
                Bs{r,d1}{k2} = repmat( k2, size(As{r,d1}{k2},1), 1 );
                
                Le{r,d1}{k2} = FEAT.SINGLE{r,k}{c}.loc( FEAT.SINGLE{r,k}{c}.cycleID<=cycleNUM_E );
            end
            Ae{r,d1} = cell2mat(Ae{r,d1}).';
            Be{r,d1} = cell2mat(Be{r,d1}).';
            As{r,d1} = cell2mat(As{r,d1}).';
            Bs{r,d1} = cell2mat(Bs{r,d1}).';
            
            Le{r,d1} = cell2mat(Le{r,d1}).';

            % split the data into folds
            for k = 1:authFoldN
                foldEleIdxB{k,r,d1} = ismember( Be{r,d1}, cohortIdx{k} );
                foldExeIdxB{k,r,d1} = ~foldEleIdxB{k,r,d1};
                foldEleIdx{k,r,d1}    = find( foldEleIdxB{k,r,d1} ).';
                foldEleNum(k,r,d1)    = length( foldEleIdx{k,r,d1} );
            end
            % for weight computation
            [~,spNum_s_r] = groupIndex( Bs{r,d1}, person_total );
            spNum_s(c,r,:,d1) = reshape( spNum_s_r, [ 1 1 person_total ] );
        end
    end
    
    % normalize before clustering
    for d1 = 1:2
        for r=1:record_total
            Ae{r,d1} = bsxfun( @rdivide, Ae{r,d1}, sqrt( sum( Ae{r,d1}.*Ae{r,d1}, 1 ) ) );
        end
    end
    
    % cluster the data
    % -- initialize the kmean seeds
    rng(randSeed);
    initCenterIdent = cell(record_total,2); % dictionary
    initCenterAuth  = cell(authFoldN, record_total,2); % dictionary
    numWordsIdent   = zeros(record_total,2);
    numWordsAuth    = zeros(authFoldN, record_total,2);
    for d1 = 1:2
        for r=1:record_total
            % for identification
            numWordsIdent(r,d1)   = getWordNum( wordNumBase, person_total );
            numWordsIdent(r,d1)   = max( ceil( numWordsIdent(r,d1)*portionSP ), 10 );
            initCenterIdent{r,d1} = randperm( size(Ae{r,d1},2), numWordsIdent(r,d1) );
            % for authentication
            for k = 1:authFoldN
                numWordsAuth(k,r,d1)   = getWordNum( wordNumBase, authFoldEleN(k) );
                initCenterAuth{k,r,d1} = randperm( foldEleNum(k,r,d1), numWordsAuth(k,r,d1) );
            end
        end
    end
    % -- perform kmeans, generate words
    WordsIdent = cell(record_total,2);
    ClusterLabelIdent = cell(record_total,2);
    
%     for r=1:record_total
    parfor r=1:record_total
        WordsIdent_r = cell(1,2);
        ClusterLabelIdent_r = cell(1,2);
        for d1 = 1:2
            % for identification
            [ClusterLabelIdent_r{1,d1}, WordsIdent_r{1,d1}] = ...
                kmeans( Ae{r,d1}.', numWordsIdent(r,d1), ...
                'start', Ae{r,d1}(:,initCenterIdent{r,d1}).', ...
                'emptyaction', 'singleton', ...
                'options', statset('MaxIter',1000) );
        end
        WordsIdent_r = arrayfun( @(a) {transpose(a{1})}, WordsIdent_r );
        ClusterLabelIdent_r = arrayfun( @(a) {transpose(a{1})}, ClusterLabelIdent_r );
        WordsIdent(r,:)  = WordsIdent_r;
        ClusterLabelIdent(r,:)  = ClusterLabelIdent_r;
    end
    unWordsIdent = WordsIdent;
    
    WordsAuth  = cell(authFoldN, record_total,2);
    ClusterLabelAuth  = cell(authFoldN, record_total,2);
    
    if hasAuth
%         for r=1:record_total
        parfor r=1:record_total
            WordsAuth_r  = cell(authFoldN,1,2);
            ClusterLabelAuth_r  = cell(authFoldN,1,2);
            for d1 = 1:2

                for k = 1:authFoldN
                    % for authentication
                    [ClusterLabelAuth_r{k,1,d1}, WordsAuth_r{k,1,d1}] = ...
                        kmeans( Ae{r,d1}(:,foldEleIdxB{k,r,d1}).', numWordsAuth(k,r,d1), ...
                        'start', Ae{r,d1}(:,foldEleIdx{k,r,d1}(initCenterAuth{k,r,d1})).', ...
                        'emptyaction', 'singleton', ...
                        'options', statset('MaxIter',1000) );
                end

            end
            WordsAuth_r  = arrayfun( @(a) {transpose(a{1})}, WordsAuth_r  );
            ClusterLabelAuth_r  = arrayfun( @(a) {transpose(a{1})}, ClusterLabelAuth_r );
            WordsAuth(:,r,:) = WordsAuth_r;
            ClusterLabelAuth(:,r,:) = ClusterLabelAuth_r;
        end
        unWordsAuth = WordsAuth;
    end
        
    
    if ~isempty( detailedOutputPath )
        save( fullfile( detailedOutputPath, ['BaseInfo-c' int2str(c) '.mat'] ), ...
            'Ae', 'Be', 'As', 'Bs', 'Le', ...
            'WordsIdent', 'ClusterLabelIdent' ...  
            ... % , 'WordsAuth', 'ClusterLabelAuth'
            );
    end
    
    fprintf(1,'+');
    
    % normalize data

    for r=1:record_total
        for d1 = 1:2
            WordsIdent{r,d1} = ...
                bsxfun( @rdivide, WordsIdent{r,d1}, ...
                sqrt( sum( WordsIdent{r,d1}.*WordsIdent{r,d1}, 1 ) ) );
            if hasAuth
                for k = 1:authFoldN
                    WordsAuth{k,r,d1} = ...
                        bsxfun( @rdivide, WordsAuth{k,r,d1}, ...
                        sqrt( sum( WordsAuth{k,r,d1}.*WordsAuth{k,r,d1}, 1 ) ) );
                end
            end
        end
    end
    
    
    % recognize

    for er=listER  % 1:record_total
        for d1 = 1:2    % d1 for exemplars, d2 for samples
            d2 = 3-d1;
            if hasIdent
                cAeident = data2cellByIndex( Ae{er,d1},ClusterLabelIdent{er,d1},2, numWordsIdent(er,d1) );
                cBeident = data2cellByIndex( Be{er,d1},ClusterLabelIdent{er,d1},2, numWordsIdent(er,d1) );
                WordsIdent_er_d   = WordsIdent{er,d1};
                unWordsIdent_er_d = unWordsIdent{er,d1};
            end
            if hasAuth
                cAeauth = cell( authFoldN, 1 );
                cBeauth = cell( authFoldN, 1 );
                dClusterLabelExe = cell( authFoldN, 1 );
                
                for k = 1:authFoldN
                    cAeauth{k} = data2cellByIndex( Ae{er,d1}(:,foldEleIdxB{k,er,d1}), ...
                        ClusterLabelAuth{k,er,d1},2, numWordsAuth(k,er,d1) );
                    cBeauth{k} = data2cellByIndex( Be{er,d1}(:,foldEleIdxB{k,er,d1}), ...
                        ClusterLabelAuth{k,er,d1},2, numWordsAuth(k,er,d1) );
                    
                    ClusterLabelExe_k   = knnclassify( Ae{er,d1}(:,foldExeIdxB{k,er,d1}).', ...
                        unWordsAuth{k,er,d1}.', 1:numWordsAuth(k,er,d1), 1 );
                    dClusterLabelExe{k} = data2cellByIndex(ClusterLabelExe_k.', ...
                        Be{er,d1}(:,foldExeIdxB{k,er,d1}), 2, person_total_both);
                end
                dAeer = data2cellByIndex(Ae{er,d1}, Be{er,d1},2, person_total);
            end
            for sr=1:record_total 
                % For both Identification & Authentication
                dAssr = data2cellByIndex(As{sr,d2},Bs{sr,d2},2,person_total);
                % Identification
                if hasIdent
                    Rpar = zeros(person_total_both,person_total);
                    parfor sk=1:person_total_both
%                     for sk=1:person_total_both
                        clusterWeight = ones( size(WordsIdent_er_d,2), size(dAssr{sk},2) );
                        if chosenClusterNum<Inf
                            switch srcType
                                case 1
                                    clusterWeight = omp( WordsIdent_er_d, dAssr{sk}, ...
                                        [], chosenClusterNum );
                                case 2
                                    for pk=1:size(dAssr{sk},2)
                                        clusterWeight(:,pk) = omp2( WordsIdent_er_d, dAssr{sk}(:,pk), ...
                                            [], clusterCodingTol * norm(dAssr{sk}(:,pk)) );
                                    end
                                case 3
                                    for pk=1:size(dAssr{sk},2)
                                        awDelta = bsxfun( @minus, unWordsIdent_er_d,dAssr{sk}(:,pk) );
                                        [~,bestIdx] = sort( sum( awDelta.*awDelta, 1 ) );
                                        clusterWeight(bestIdx(1:chosenClusterNum),pk) = 1;
                                    end
                                otherwise
                                    error( 'Unrecognized srcType' );
                            end
                        end
                        chosenClusterIdxB = full( clusterWeight ~=0 );
                        Rpar(sk,:) = zgSPScore2( dAssr{sk}, chosenClusterIdxB, ...
                            cAeident, cBeident, person_total, tolParameter, srcType, softmaxSigma2, ...
                            detailedOutputPath, sprintf('ident-c%d-d%d-sr%d-er%d-sk%03d.mat',c,d1,sr,er,sk) );
                    end
                    R.IDENT( c, sr, er, :, :, d1 ) = shiftdim( Rpar, -3 );
                end

                % authentication
                if hasAuth
                    for k=1:authFoldN
                        WordsAuth_k_er_d     = WordsAuth{k,er,d1};
                        cAauth_k             = cAeauth{k};
                        cBauth_k             = cBeauth{k};
                        
                        numWordsAuth_k_er_d1 = numWordsAuth(k,er,d1);
                        authFoldSamN_k       = authFoldSamN(k);
                        samIdx_k             = samIdx{k};
                        exeIdx_k             = exeIdx{k};
                        cohortIdx_k          = cohortIdx{k};
                        
                        dClusterLabelExe_k1   = dClusterLabelExe{k}(exeIdx_k);
                        dAeer1 = dAeer(exeIdx_k);
                        
                        Rpar = zeros( length(cohortIdx_k)+1, authFoldSamN_k , length(exeIdx_k) );
                        
                        parfor ek0=1:length(exeIdx_k)
%                         for ek0=1:length(exeIdx_k)
                            ek = exeIdx_k(ek0);
                            % integrate exemplar and cohort
                            singleExemplarCluster = data2cellByIndex( dAeer1{ek0}, ...
                                dClusterLabelExe_k1{ek0}, 2, numWordsAuth_k_er_d1 );
                            singleExemplarClusterLen = cellfun( @(x) size(x,2), singleExemplarCluster );
                            cAauth_e = arrayfun( @(a1,a2) {[a1{1},a2{1}]}, cAauth_k, singleExemplarCluster );
                            cBauth_e = arrayfun( @(b1,l1) {[b1{1},repmat(ek,1,l1)]}, cBauth_k, singleExemplarClusterLen );
                            % authenticate
%                             dAsr1 = dAsr(samIdx_k);

                            dAssrC = dAssr;
                            samIdx_kC = samIdx_k;
                            for sk0=1:authFoldSamN_k
                            	sk = samIdx_kC(sk0);
                                clusterWeight = ones( size(WordsAuth_k_er_d,2), size(dAssrC{sk},2) );
                                if chosenClusterNum < Inf
                                    switch srcType
                                        case 1
                                            clusterWeight = omp( WordsAuth_k_er_d, dAssrC{sk}, ...
                                                [], chosenClusterNum );
                                        otherwise
                                            error( 'Unsupported srcType for authentication.' );
                                    end
                                end
                                chosenClusterIdxB = full( clusterWeight ~=0 );
                                R_sk = zgSPScore2( dAssrC{sk}, chosenClusterIdxB, ...
                                    cAauth_e, cBauth_e, person_total_both, chosenExemplarNum, srcType, softmaxSigma2 );
                                Rpar(:,sk0,ek0) = R_sk([ek cohortIdx_k]);

                            end
                        end
                        R.AUTH{k}( c, sr, er, :, :, :, d1 ) = shiftdim( Rpar, -3 );
                    end
                end

            end
        end
        fprintf(1,'.');
    end
end


fprintf(1,'\n');

end

%%
function R = zgSPScore2( X, chosenClusterIdxB, cA, cB, N, param1, srcType, softmaxSigma2, detailedOutputPath, detailedOutputFileName )


switch srcType
    case 1
        chosenExemplarNum = param1;
    case 2
        exemplarCodingTol = param1;
    case 3
        chosenExemplarNum = param1;
    otherwise
        error( 'Unrecognized srcType' );
end

spN   = size(X,2);

P = zeros(N,spN);

GaNorm1 = zeros(N,spN);

for k = 1:spN
    
    % coding
    cAk = cell2mat( cA( :, chosenClusterIdxB(:,k) ) );
    cBk = cell2mat( cB( :, chosenClusterIdxB(:,k) ) );
    switch srcType
        case 1
            ga = omp( cAk, X(:,k), [], chosenExemplarNum );
        case 2
            ga = omp2( cAk, X(:,k), [], exemplarCodingTol * norm(X(:,k)) );
        case 3
            ga = omp( cAk, X(:,k), [], chosenExemplarNum );
        otherwise
            error( 'Unrecognized srcType' );
    end
    ga = full(ga);
    
    % src
    cAk_parted = data2cellByIndex( cAk, cBk, 2, N );
    ga_parted  = data2cellByIndex( ga , cBk, 1, N );
    
    GaNorm1(:,k) = cellfun( @(x) sum(abs(x)), ga_parted ); % for detailed output ONLY
    
    H = arrayfun( @(x, A) { A{1}*x{1} }, ga_parted, cAk_parted.' );
    % ---- not necessary
%     isHempty = cellfun( @isempty, H );
%     H(isHempty) = { zeros(size(X,1),1) };
    % ------------------
    P(:,k) = cellfun(  @(h) sum( abs(X(:,k)-h ) ), H );

end

useMP = 0;

if useMP
    P = -mp(P)/mp(2*softmaxSigma2);
    P = P - repmat( max(P,[],1), size(P,1), 1 );
    P = exp(P);
    Rarr = P .* repmat( 1./sum(P,1), size(P,1), 1 );
else
    P = -P/(2*softmaxSigma2);
    P = bsxfun( @minus, P, max(P,[],1) );
    P = exp(P);
    Rarr = bsxfun( @times, P, 1./sum(P,1) );
end

if exist( 'detailedOutputPath', 'var' )
    if ~isempty( detailedOutputPath )
        save( fullfile( detailedOutputPath, detailedOutputFileName ) , 'Rarr', 'chosenClusterIdxB', 'GaNorm1' );
    end
end

Rarr = log(Rarr);
R_unnorm = sum(Rarr,2);
% R_shift  = R_unnorm - min(R_unnorm);
% R_shift_exp = exp(R_shift);
% R = log( R_shift_exp/sum(R_shift_exp) );
% W = repmat( sum( exp( R_unnorm ) ), size(R_unnorm) );
R = double( R_unnorm );

end

%%
function dA = data2cellByIndex(A,B,dim, varargin)

if ~exist('dim','var')
    dim = 1;
end

[ sortedIdx, cellLen ] = groupIndex(B, varargin{:} );
S = num2cell( size(A) );
S{dim} = cellLen;
I = arrayfun( @(x) {1:x} , size(A) );
I{dim} = sortedIdx;
dA = mat2cell( A(I{:}), S{:} );

end

function [ sortedIdx, cellLen ] = groupIndex(B, groupNum)

B  = reshape(B,1,numel(B));
ub = unique(B);

if ~exist( 'groupNum', 'var' )
    groupNum = max(ub);
end

[sortedL,sortedIdx] = sort(B);
deltaL   = diff( sortedL, 1 );
cellLen0 = diff( [0, find(deltaL), length(B)], 1);

cellLen     = zeros( 1, groupNum );
cellLen(ub) = cellLen0;

end

%% 
function m = getWordNum( wordNumBase, n ) % this is only a heuristic
    m = floor( wordNumBase * max( log2(n/40), 1 ) );
    m = round(m/5)*5;
end
