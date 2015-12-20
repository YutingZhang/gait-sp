function A = statistic_detailed( infoPath )

load ../results/spsrc/feat.mat
load gait
load score_param

idxChosenPerson = true( size( DATA.BOTH,2 ),1 );
idxChosenPerson( score_param.CorruptedPerson ) = false;

DATA.BOTH = DATA.BOTH( :,idxChosenPerson,: );
FEAT.BOTH = FEAT.BOTH( :,idxChosenPerson,: );

numChannel = 5;
numRecord   = size( DATA.BOTH,1 );
numBoth     = size( DATA.BOTH,2 );

% 

C = {};
for k=1:numChannel
    Cs = nchoosek( 1:numChannel, k );
    C  = [C ; mat2cell( Cs, ones(1,size(Cs,1)) , size(Cs,2) )];
end

A.numGSigCorrect = zeros(numChannel,2,numRecord,numRecord,numBoth);
A.numGSigTotal   = zeros(numChannel,2,numRecord,numRecord,numBoth);

A.meanLocationalRelevance = cell(numChannel,2,numRecord,numRecord,numBoth);

for c = 1:numChannel
    load( fullfile( infoPath, ...
        sprintf('BaseInfo-c%d.mat',c) ) );
    for d1 = 1:2
        d2 = 3-d1;
        for er = 1:numRecord
            clusterLocation = data2cellByIndex( Le{er,d1}, ClusterLabelIdent{er,d1}, 2, max( ClusterLabelIdent{er,d1} ) );
            
            for sr = 1:numRecord
                for sk = 1:numBoth
                    clear Rarr chosenClusterIdxB
                    
                    load( fullfile( infoPath, ...
                        sprintf('ident-c%d-d%d-sr%d-er%d-sk%03d.mat',c,d1,sr,er,sk) ) );
                    
                    % single G-Signature correct
                    [ ~, L ] = max( Rarr, [] , 1 );
                    A.numGSigCorrect(c,d1,sr,er,sk) = sum( L == sk );
                    A.numGSigTotal(c,d1,sr,er,sk)   = length(L);

                    % locational relevance
                    LR = zeros( 1, size(chosenClusterIdxB,2) );
                    for u = 1:size(chosenClusterIdxB,2)
                        sigLoc = FEAT.BOTH{sr,sk,d2}{c}.loc( u );
                        subLoc = cat( 2, clusterLocation{chosenClusterIdxB(:,u)} );
                        dLoc = abs(subLoc-sigLoc);
                        dLoc = min(dLoc,1-dLoc);
                        LR(u) = mean(dLoc);
                    end
                    A.meanLocationalRelevance{c,d1,sr,er,sk} = LR;
                end
            end
        end
    end
end

A.correctNum4Ch  = zeros( length( C ) );
A.totalNum4Ch    = zeros( length( C ) );
A.correctRate4Ch = zeros( length( C ) );

A.locationalHistEdges = 0:0.05:0.5;
A.locationalHist = histc( cat( 2, A.meanLocationalRelevance{:} ), 0:0.05:0.5 );

for k = 1:length( C )
    A.correctNum4Ch  = sum( sum( A.numGSigCorrect( C{k}, : ) ) );
    A.totalNum4Ch    = sum( sum( A.numGSigTotal( C{k}, : ) ) );
    A.correctRate4Ch = A.correctNum4Ch ./ A.totalNum4Ch ;
end



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
