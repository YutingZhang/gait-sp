function [score W] = zgScore_AllCycles( D1, D2, param )
% D1, D2 - two sets of features: D1 - sample, D2 - exemplar
% param: .p, .p_t, .p_v, .kernel_type, .cycleNUM

CN = length(D1);
S = zeros(1,CN);
L = zeros(1,CN);

[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

for c=1:CN
    d1 = D1{c};
    d2 = D2{c};
    n1 = min(cycleNUM_S,size(d1,2));
    n2 = min(cycleNUM_E,size(d2,2));
    dis = zeros( n1, n2 );
    for k=1:n1
        d1k = d1(:,k);
        d2c = d2(:,1:n2);
        delta = d2c - repmat( d1k, 1, n2 );
        dis(k,:) = sqrt(sum(delta.*delta));
    end
    S(c) = -min(reshape(dis,numel(dis),1));   % make it negative
    L(c) = 1;
end
score = S;
W = L./sum(L);

end

