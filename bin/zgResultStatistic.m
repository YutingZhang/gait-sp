function S = zgResultStatistic( R, W, PersonNum, param )
%
% param.   PreNormFunName: 'none', 'exp'
% param.   CorruptedPerson: ([] defaultly)
% param.   CohortFoldN:     ( 3 defaultly )

if ~exist('param','var')
    param = [];
end

% W.IDENT = ones( size(R.IDENT) );
% for k=1:length(R.AUTH_NORM)
%     W.AUTH_NORM{k} = ones( size(R.AUTH_NORM{k}) );
% end

is_auth_folded = false;
is_authNorm_ok = false;
is_authNorm_same = false;
if isstruct( R )
    Ri = R.IDENT;
    Wi = W.IDENT;
    if isfield( R, 'AUTH' )
        Ra = R.AUTH;
        Wa = W.AUTH;
        is_auth_folded = true;
    end
    if isfield( R, 'AUTH_NORM' )
        Ran = R.AUTH_NORM;
        Wan = W.AUTH_NORM;
        is_authNorm_ok = true;
    end
else
    Ri = R;
    Wi = W;
end
clear R W

bothNum   = PersonNum.both;
singleNum = PersonNum.single;
totalNum  = bothNum+singleNum;
if totalNum ~= size(Ri,5);
    error('PersonNum is not consistent with Ri');
end

% channel, sample record number, exemplar record number,
% sample persons, exemplar persons, days
rs = size(Ri);
CaseNum    = rs(2)*rs(3)*rs(6); % [sample record number, exemplar record number, days]
ChannelNum = rs(1);

Ri( :, :,:,(bothNum+1):end,(bothNum+1):end, : ) = -inf;
Wi( :, :,:,(bothNum+1):end,(bothNum+1):end, : ) = 1;

C = {};
for k=1:ChannelNum
    Cs = nchoosek( 1:ChannelNum, k );
    C  = [C ; mat2cell( Cs, ones(1,size(Cs,1)) , size(Cs,2) )];
end

if isfield( param, 'CohortFoldN' )
    foldN = param.CohortFoldN;
else
    foldN = 3;
end
[~, foldRange] = zgFoldPartition( bothNum, foldN );
cohortIDX_A = false( bothNum, foldN );
testIDXc_A = false( bothNum, foldN );
testIDXr_A = false( totalNum, foldN );
for k=1:foldN
    cohortIDX = false(bothNum,1);
    cohortIDX(foldRange(k,1):foldRange(k,2)) = true;
    testIDXc = ~cohortIDX;
    testIDXr = true(totalNum,1);
    testIDXr(1:bothNum) = and( testIDXr(1:bothNum) , testIDXc );
    
    cohortIDX_A(:,k) = cohortIDX;
    testIDXc_A(:,k)  = testIDXc;
    testIDXr_A(:,k)  = testIDXr;
end


S.IDENT = cell( length(C), 1 );
S.IDENT_PD = cell( length(C), 1 );
S.AUTH = cell( length(C), 1 );
S.AUTH_NORM = cell( length(C), 1 );

for c=1:length(C)
    C0 = C{c};
    
    RS = WeigthedResult_WithReshaping( ...
        Ri( C0, :, :, :, :, : ), Wi( C0, :, :, :, :, : ) );
    
    % identification
    [~, IDX] = sort(RS(1:bothNum,:,:),2, 'descend');  % exemplar
    RANK_IND = (IDX == repmat((1:bothNum).',[1 totalNum CaseNum]));
    T = sum( RANK_IND,3 );  % correct rate at different rank (dim 2 for rank)
    S.IDENT{c} = sum( T, 1 );
    S.IDENT_PD{c} = sum(repmat(1:totalNum,[bothNum 1]).*T,2);
    
    % non-normalized authentication
    if is_auth_folded
        RSa = WeigthedResult_WithReshaping4Fold( ...
                Ra, Wa, C0 );
        S.AUTH{c} = auth_result( RSa );
        if length( size( Ra{1} ) ) == 7 % for MAP fusion. AUTH_NORM == AUTH
%             is_authNorm_ok = true;
            is_authNorm_same = true;
%             Ran = Ra; Wan = Wa;
        end
    else
        RSa = RS( :, 1:bothNum, : );
        S.AUTH{c} = auth_result( RSa );
    end
    
    % normalized authentication
    
    if is_authNorm_same
        
        S.AUTH_NORM{c} = S.AUTH{c};
        
    else
    
        if is_authNorm_ok
            RSaC = WeigthedResult_WithReshaping4Fold( ...
                Ran, Wan, C0 );
            S.AUTH_NORM{c} = auth_result( RSaC );
        else

            if size(RSa,1)<totalNum     % in case that RS ONLY contains results for IDENT instead of the FULL results
                S.AUTH_NORM = S.AUTH;
            else

                if ~isfield( param, 'PreNormFunName' )
                    param.PreNormFunName = 'none';
                end
                [ NormCallBack RSaN ] = zgScoreNormInit( RSa , param.PreNormFunName );

                RSaC = cell(foldN,1);
                for k=1:foldN
                    cohortIDX = cohortIDX_A(:,k);
                    testIDXc  = testIDXc_A(:,k);
                    testIDXr  = testIDXr_A(:,k);
                    RSaN_t = RSaN(testIDXr,testIDXc,:);
                    RSaN_c = RSaN(testIDXr,cohortIDX,:);
                    RSaC{k}  = feval( NormCallBack, ...
                        log(RSaN_t)-log(repmat( mean(RSaN_c,2), [ 1 size(RSaN_t,2) 1 ] )) ...
                        );
                end
                S.AUTH_NORM{c} = auth_result( RSaC );
            end

        end
        
    end
    
    
    
end

S.AUTH      = cell2mat(S.AUTH);
S.AUTH_NORM = cell2mat(S.AUTH_NORM);

end

function RSa = WeigthedResult_WithReshaping4Fold( Ra, Wa, C0 )
    foldN = length(Ra);
    RSa = cell( foldN, 1 );
    for k=1:foldN
        fts = size(Ra{k});
        RSa{k} = WeigthedResult_WithReshaping( ...
            reshape( Ra{k}( C0, : ), [length(C0),fts(2:end)]), ...
            reshape( Wa{k}( C0, : ), [length(C0),fts(2:end)]) );
    end
end

function RS = WeigthedResult_WithReshaping( R0, W0 )



if length( size(R0) ) == 7  % MAP fusion
    useMP = 0;
    
%     R0 = R0./50;
    
    if useMP
        % advanpix doesn't support multi-dim array ...
        RSall    = sum(R0.*W0,1);
        origSize = size( RSall );
        RSall    = permute( RSall, [4 1 2 3 5 6 7] );
        RSall    = RSall(:,:);
        
        subTaskScale = 1e4;
        batchNum = floor( size( RSall, 2) / subTaskScale );
        batchMod = mod( size( RSall, 2), subTaskScale );
        cellSplitter = repmat( subTaskScale, 1, batchNum );
        if batchMod
            cellSplitter = [ cellSplitter, batchMod ];
        end
        
        RSall = mat2cell( RSall, size(RSall,1), cellSplitter );
        RS = cell( size(RSall) );
        
        parfor k = 1:numel(RSall)
            RSall_k = mp( RSall{ k } );
            maxRSall_k = max(RSall_k,[],1);
            for i = 1:size( RSall_k, 1 )
                RSall_k(i,:) = RSall_k(i,:) - maxRSall_k;
            end
            RS_k    = RSall_k(1,:) - log( sum( exp(RSall_k), 1 ) );
            RS{k}   = double( RS_k );
        end
        clear RSall
        RS = cell2mat( RS );
        
        RS    = reshape( RS, [1 origSize(1:3) origSize(5:7)] );
        RS    = permute( RS, [ 2 3 4 1 5 6 7 ] );
        RS    = permute( RS, [ 5 6 2 3 7 1 4 ]);
        RS    = RS(:,:,:);
        

    else
        RSall = sum(R0.*W0,1);
        % R: channel, sr, er, cohortIdx, sk, ek, d1
        RSall = bsxfun( @minus, RSall, max( RSall, [] , 4 ) ) + 600; % 600 is constant no more than 700 (the max)

        RS = RSall(:,:,:,1,:,:,:) - log( sum( exp(RSall), 4 ) );


    %     RS = permute(RS, [ 1 2 3 5 6 7 4 ]);
    %     RS = permute(RS, [ 4 5 2 3 6 1 ]);
        RS = permute(RS, [ 5 6 2 3 7 1 4 ]);
        RS = RS(:,:,:);
    end
    fprintf( 1, '.' );
else  % normal case
%     rs = size(R0);
%     CaseNum = rs(2)*rs(3)*rs(6);
    RS = sum(R0.*W0,1); %./sum(W0,1);
    RS = permute(RS, [ 4 5 2 3 6 1 ]);
    % sample number, exemplar number, case number
%     RS = reshape( RS, [rs(4), rs(5), CaseNum] );
    RS = RS(:,:,:);
end

end

function AUTH = auth_result( RSa )

if ~iscell( RSa )
    RSa = {RSa};
end

INDICATOR = cell( length( RSa ) , 1 );
Q = cell( length( RSa ) , 1 );

for k=1:length( RSa )
    RSa0 = RSa{k};
    rn = size( RSa0, 1 );
    cn = size( RSa0, 2 );
    FIT_IDX0 = ((1:cn)-1)*(rn+1)+1;
    FIT_IDX  = false(rn*cn,1); FIT_IDX(FIT_IDX0)  = true;
    Q{k} = reshape( RSa0, numel(RSa0), 1 );
    INDICATOR{k}  = repmat( FIT_IDX, size(RSa0,3), 1 );
end

Q = cell2mat( Q );
INDICATOR = cell2mat( INDICATOR );

[ AUTH.SORTED_SCORE SIDX ] = sort( Q, 'ascend' );
AUTH.FIT_INDICATOR = INDICATOR(SIDX);

AUTH.FALSE_REJECT = cumsum(AUTH.FIT_INDICATOR);
AUTH.FALSE_ACCEPT = cumsum(~AUTH.FIT_INDICATOR(end:-1:1));
AUTH.FALSE_ACCEPT = AUTH.FALSE_ACCEPT(end:-1:1);

AUTH.ACCEPT_TOTAL = sum(INDICATOR);
AUTH.REJECT_TOTAL = length(AUTH.SORTED_SCORE)-AUTH.ACCEPT_TOTAL;

pFR = AUTH.FALSE_REJECT/AUTH.ACCEPT_TOTAL;
pFA = AUTH.FALSE_ACCEPT/AUTH.REJECT_TOTAL;
eerIDX = find(pFR>pFA,1);

% get EER from interpolation
Dy2 = pFR(eerIDX)-pFR(eerIDX-1);
Dy1 = pFA(eerIDX)-pFA(eerIDX-1);
AUTH.EER  = ( pFR(eerIDX)*Dy1 - pFA(eerIDX)*Dy2 ) / (Dy1-Dy2);

% sampling

if 1
    [FR1 FA1] = simplifyPlot( AUTH.FALSE_REJECT, AUTH.FALSE_ACCEPT, 1000 );
    [FA2 FR2] = simplifyPlot( AUTH.FALSE_ACCEPT(end:-1:1), AUTH.FALSE_REJECT(end:-1:1), 1000 );
    FRx = [FR1;FR2]; FAx = [FA1;FA2];
    [ AUTH.FALSE_REJECT idx2 ] = sort( FRx,'ascend' );
    AUTH.FALSE_ACCEPT = FAx( idx2 );
end

% clear useless field
if 1
    AUTH = rmfield( AUTH, 'FIT_INDICATOR' );
    AUTH = rmfield( AUTH, 'SORTED_SCORE' );
end

end

function [ Xi Yi ] = simplifyPlot( X, Y, XsliceNum )

% make X distinct
[ X_S, LAST ]     = unique(X,'last');
Y_S = ( Y( [1; LAST(1:end-1)+1] ) + Y(LAST) )/2;
if X(1)~=X_S(1)
    X = [ X(1) ; X_S ];
    Y = [ Y(1) ; Y_S ];
end
if X(end)~=X_S(end)
    X_S = [ X_S; X(end) ];
    Y_S = [ Y_S; Y(end) ];
end

% interplate
x1 = min(X_S); x2 = max(X_S);
Xi = x1:((x2-x1)/XsliceNum):x2;
if size(X,1)>1
    Xi = Xi.';
end
Yi = interp1( X_S, Y_S, Xi );

end

