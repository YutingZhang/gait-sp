function [score W] = zgScore_GSig( D1, D2, param )
% D1, D2 - two sets of features: D1 - sample, D2 - exemplar
% param: .p, .p_t, .p_v, .kernel_type, .cycleNUM

CN = length(D1);
S = zeros(1,CN);
L = zeros(1,CN);

if ~isfield(param,'nearest_k')
    param.nearest_k = inf;
end

[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

if (    ( cycleNUM_E<inf && ~isfield( D2{1}, 'cycleID' ) ) || ...
        ( cycleNUM_S<inf && ~isfield( D1{1}, 'cycleID' ) ) )
    error('cycleID is not available, but needed.');
end

for c=1:CN
    if cycleNUM_S<inf
        A1 = sampleFirstCycles(D1{c}, cycleNUM_S);
    else
        A1 = D1{c};
    end
    if cycleNUM_E<inf
        A2 = sampleFirstCycles(D2{c}, cycleNUM_E);
    else
        A2 = D2{c};
    end
    S(c) = zgGSigMatch( A1.descriptor, A1.loc, A1.scale, ...
        A2.descriptor, A2.loc, A2.scale, ...
        param.p, param.p_t, param.p_v, param.kernel_type, param.nearest_k );
    L(c) = length(D1{c}.loc);
end
score = S;
W = L./sum(L);
% score = sum(S.*L)/sum(L);
% score = exp(score);

end

function B = sampleFirstCycles( A, cycleNUM )

IDX = ( A.cycleID <= cycleNUM );
B.descriptor = A.descriptor(IDX,:);
B.loc   = A.loc(IDX);
B.scale = A.scale(IDX);

end
