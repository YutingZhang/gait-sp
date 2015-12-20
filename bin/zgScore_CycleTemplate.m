function [score W] = zgScore_CycleTemplate( D1, D2, param )
% D1, D2 - two sets of features: D1 - sample, D2 - exemplar
% param: .p, .p_t, .p_v, .kernel_type, .cycleNUM


l_order = 2;
if exist('param','var')
    if isfield( param, 'l_order' )
        l_order = param.l_order;
    end
end

CN = length(D1);
S = zeros(1,CN);
L = zeros(1,CN);

[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

% compute features
F1 = zgTemplateFeatures(D1,cycleNUM_S);
F2 = zgTemplateFeatures(D2,cycleNUM_E);

% compare
for c=1:CN
    S(c) = - sum(abs(F2{c}-F1{c}).^(l_order)).^(1/l_order); % make it negative
    L(c) = 1;
end
score = S;
W = L./sum(L);

end

function D = zgTemplateFeatures(CY, cycleNUM)

CN = length( CY );
D=cell(size(CY));

for c=1:CN
    CS = CY{c};
    n = min(size(CS,2),cycleNUM);
    A = mean( CS(:,1:n), 2 );
    D{c} = A;
end

end

