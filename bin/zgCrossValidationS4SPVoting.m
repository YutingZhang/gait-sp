function [R W paramPersonNum] = zgCrossValidationS4SPVoting( FEAT , score_fun , param, options )
% R = zgCrossValidationS( FEAT , score_fun , param )
if ischar(FEAT)
    load( FEAT );
end

paramPersonNum.both  = size(FEAT.BOTH,2);
paramPersonNum.single = size(FEAT.SINGLE,2);

%% Override score_fun
% score_fun = @zgScoreArr_SPVoting;

%% Init num
person_total_both = size(FEAT.BOTH,2);
record_total = size(FEAT.BOTH,1);
person_total_single = size(FEAT.SINGLE,2);
person_total = person_total_both+person_total_single;

%% Parsing options
authFoldN = 3;
if exist( 'options', 'var' )
    if isfield( options, 'authFoldN' )
        authFoldN = options.authFoldN;
    end
end
[ authFoldEleN, authFoldRange ] = zgFoldPartition( person_total_both, authFoldN );

%% Parsing Param
if isfield( param, 'norm_order' )
    normOrder = param.norm_order;
    fprintf(1,'norm_order = %f\n', normOrder);
else
    normOrder = 1;
end

hasAuth = 1;
if isfield( param, 'hasAuth' )
    hasAuth = param.hasAuth;
end


%% Process

% channel, sample record number, exemplar record number,
% sample persons, exemplar persons, days
R.IDENT = zeros( [ 5, record_total , record_total , ...
    person_total_both , person_total, 2 ] );
ts = size(R.IDENT);
channelN = ts(1);
W.IDENT  = zeros(ts);

IR1 = zeros( channelN, 1, 1, person_total_both, person_total );
IW1 = zeros( channelN, 1, 1, person_total_both, person_total );

R.AUTH   = cell(authFoldN,1); W.AUTH   = cell(authFoldN,1);
R.AUTH_NORM = cell(authFoldN,1); W.AUTH_NORM = cell(authFoldN,1);
cohortIDXArr = false(person_total,authFoldN);
authFoldSamN = zeros(authFoldN,1);
authFoldExeN = zeros(authFoldN,1);
cohortIDXNum = cell(authFoldN,1);
samIDXNum = cell(authFoldN,1);
for k = 1:authFoldN
    authFoldSamN(k) = person_total-authFoldEleN(k);
    authFoldExeN(k) = person_total_both-authFoldEleN(k);
    R.AUTH{k} = zeros( [ 5, record_total , record_total , ...
        authFoldSamN(k) , person_total_both-authFoldEleN(k), 2 ] );
    fts = size(R.AUTH{k});
    W.AUTH{k} = zeros(fts);
    R.AUTH_NORM{k} = zeros(fts);
    W.AUTH_NORM{k} = zeros(fts);
    
    cohortIDXNum{k} = authFoldRange(k,1):authFoldRange(k,2);
    cohortIDXArr( cohortIDXNum{k}, k ) = true;
    samIDXNum{k} = find( ~cohortIDXArr( :, k ) ).';
    
end
exeIDXArr = ~cohortIDXArr(1:person_total_both,:);

AR1A_eN = max( authFoldExeN(k) );
AR1A = zeros( channelN, 1, 1, person_total, person_total_both, authFoldN );
AR1A_norm = zeros( size(AR1A) );
AR1A_nw = zeros( size(AR1A) );


for d1 = 1:2
    d2 = 3-d1;
    for er=1:record_total
        %% Identification & authentication
        EXE = [ FEAT.BOTH(er,:,d1) , FEAT.SINGLE(er,:)];
        for sr=1:record_total
            SAM = [ FEAT.BOTH(sr,:,d2) , FEAT.SINGLE(er,:) ];

            parfor sk=1:person_total
%             for sk=1:person_total
                sample = SAM{sk};
                if sk<=person_total_both
                    % identification
                    DIST = zgSPCompare( sample, EXE, param );
                    [R0, W0] = zgSPVotingScore( DIST, 1:person_total, param);
                    IR1(:,1,1,sk, :) = reshape( R0, [ size(R0,1), 1,1,1, size(R0,2) ]);
                    IW1(:,1,1,sk, :) = reshape( W0, [ size(R0,1), 1,1,1, size(W0,2) ]);
                else
                    % authentication
                    DIST = zgSPCompare( sample, FEAT.BOTH(er,:,d1), param );
                end
                
                if hasAuth
                    % compute weight
                    ACnw0 = zgSPChannelWeight( sample );
                    % authentication
                    authFoldEleN_p = authFoldEleN;
                    samIDXNum_p    = samIDXNum;
                    cohortIDXNum_p = cohortIDXNum; 
                    authFoldExeN_p = authFoldExeN;
                    exeIDXArr_p = exeIDXArr;
                    AR1F = zeros( channelN, 1, 1, 1, AR1A_eN, authFoldN );
                    AR1F_norm = zeros( size(AR1F) );
                    AR1F_nw   = zeros( size(AR1F) );
                    for k = 1:authFoldN
                        if ismember(sk, cohortIDXNum_p{k})
                            continue;
                        end
                        eN = authFoldExeN_p(k);
                        AC = zeros( channelN, authFoldEleN_p(k)+1, eN );
                        eIDX = samIDXNum_p{k};
                        for eko=1:eN
                            ekt = eIDX(eko);
                            [R0,~] = zgSPVotingScore( DIST, [ekt cohortIDXNum_p{k}], param );
                            AC(:,:,eko) = R0;
                        end
                        ACc    = AC( :,1,: );
                        ACcl   = ACc.^normOrder; % ACc is positive, no need for abs
                        ACs    = sum(ACcl,2);
                        ACs = ACs + eps .*(ACs==0);
                        ACnorm = (ACcl./ACs).^(1/normOrder);
    %                     ACnw   = ACc./repmat(sum(ACc,1),[channelN, 1, 1]);
                        ACnw   = repmat( ACnw0, [1, 1, eN] );
                        AR1F( :, 1, 1, 1, exeIDXArr_p(:,k), k )      = reshape( ACc,    [channelN 1 1 1 eN] );
                        AR1F_norm( :, 1, 1, 1, exeIDXArr_p(:,k), k ) = reshape( ACnorm, [channelN 1 1 1 eN] );
                        AR1F_nw( :, 1, 1, 1, exeIDXArr_p(:,k), k )   = reshape( ACnw,   [channelN 1 1 1 eN] );
                    end
                    AR1A(:,1,1,sk, :,: ) = AR1F;
                    AR1A_norm(:,1,1,sk, :,: ) = AR1F_norm;
                    AR1A_nw(:,1,1,sk, :,: )   = AR1F_nw;
                end
            end
            
            R.IDENT(:,sr,er,:, :, d1) = IR1;
            W.IDENT(:,sr,er,:, :, d1) = IW1;
            
            if hasAuth
                for k = 1:authFoldN
                    ARC = AR1A( :, :,:, ~cohortIDXArr(:,k), exeIDXArr(:,k), k );
                    R.AUTH{k}(:,sr,er,:, :, d1) = ARC;
                    W.AUTH{k}(:,sr,er,:, :, d1) = ones(size(ARC));

                    R.AUTH_NORM{k}(:,sr,er,:, :, d1) = AR1A_norm( :, :,:, ~cohortIDXArr(:,k), exeIDXArr(:,k), k );
                    W.AUTH_NORM{k}(:,sr,er,:, :, d1) = AR1A_nw( :, :,:, ~cohortIDXArr(:,k), exeIDXArr(:,k), k );
                end
            end
            
        end
        fprintf(1,'.');
    end
end
fprintf(1,'\n');

end

%%
function sW = zgSPChannelWeight( sample )

CN = length( sample );
sW = zeros( CN,1 );
for c=1:CN
    sW(c) = size(sample{c}.descriptor,1);
end
sW = sW./sum(sW);

end

%%
function DIST = zgSPCompare( sample, exemplars, param )

if isfield( param, 'p_t' )
    p_t = param.p_t;
else
    p_t = 0.15;
end

% it's not efficient to put it here. But for extensibility
[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

CN = length(sample);
eN = length(exemplars);

DIST = cell(CN,1);

for c=1:CN
    eD = cell(eN,1);
    eI = cell(eN,1);
    eL = cell(eN,1);
    eID= cell(eN,1);
    for k=1:eN
        eD{k} = exemplars{k}{c}.descriptor;
        eL{k} = exemplars{k}{c}.loc;
        eID{k}= exemplars{k}{c}.cycleID;
        eI{k} = repmat( k, size(eD{k},1), 1 );
    end
    eD = cell2mat(eD);
    eI = cell2mat(eI);
    eL = cell2mat(eL);
    eID = cell2mat(eID);
    eaIDX = (eID<=cycleNUM_E);
    eD=eD(eaIDX,:); eI=eI(eaIDX); eL=eL(eaIDX);
    
    saIDX = ( sample{c}.cycleID<=cycleNUM_S );
    sD = sample{c}.descriptor(saIDX,:);
    sL = sample{c}.loc(saIDX);
    sDN = size(sD,1);
    DIST{c} = cell(sDN,1);
    for k=1:sDN
        dL = abs(sL(k)-eL);
        dL = min(1-dL,dL);
        aIDX = (dL<p_t);
        rN = sum(aIDX);
        pd.categoryID = eI(aIDX);
        dD = repmat(sD(k,:),rN,1)-eD(aIDX,:);
        pd.distance = sum(dD.*dD,2); % L2_norm^2
        DIST{c}{k} = pd;
    end
end

end

%%
function [ R W ] = zgSPVotingScore( DIST, exemplarIDX, param )

if isfield( param, 'k' )
    knn_k = param.k;
else
    knn_k = 1;
end

CN = length(DIST);
eN = length(exemplarIDX);

idxMAP = zeros(max(exemplarIDX),1);
idxMAP(exemplarIDX) = 1:eN;

R = zeros( CN,eN );
W = ones ( CN,eN );
for c=1:CN
    sDN = length(DIST{c});
    sI = zeros(sDN,1);
    for k=1:sDN
        pd = DIST{c}{k};
        aeIDX = ismember( pd.categoryID, exemplarIDX );
        sidx = knnFromDistance( pd.distance(aeIDX), pd.categoryID(aeIDX), knn_k );
        sI(k) = idxMAP(sidx);
    end
    R(c,:) = histc(sI,1:eN);
end

end

function testL = knnFromDistance( dist, tranL, k )

[ ~, midx ] = mink( dist, k );
sL = tranL(midx);
[uL, ~, mL] = unique(sL);
H = histc( mL , 1:length(mL) );
[maxV, maxID]= max(H);

candidateL = mL( H == maxV );
if numel(candidateL)>1
    idx = find(ismember( mL, candidateL ),1);
    maxID = mL(idx);
end

testL = uL(maxID);

end
