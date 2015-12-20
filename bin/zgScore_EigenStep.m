function [score W] = zgScore_EigenStep( D1, D2, param )
% D1, D2 - two sets of features: D1 - sample, D2 - exemplar
% param: .p, .p_t, .p_v, .kernel_type, .cycleNUM

CN = length(D1);
S = zeros(1,CN);

[ cycleNUM_S cycleNUM_E ] = zgReadCycleNumFromParam( param );

for c=1:CN
    d1 = D1{c}.Y{ min( cycleNUM_S,length(D1{c}.Y) ) };
    d2 = D2{c}.Y{ min( cycleNUM_E,length(D2{c}.Y) ) };
    switch param.cvsRecType
        case 'ident'
            y1 = d1.ident{ param.cvsDayID }.y;
            y2 = d2.ident{ param.cvsDayID }.y;
            eng = d1.ident{ param.cvsDayID }.eng;
        case 'auth'
            y1 = d1.auth{ param.cvsFoldID, param.cvsDayID }.y;
            y2 = d2.auth{ param.cvsFoldID, param.cvsDayID }.y;
            eng = d1.auth{ param.cvsFoldID, param.cvsDayID }.eng;
        otherwise
            error( 'Unknown cvsRecType' );
    end
%     dy = abs(y2-y1);
%     S(c) = sum(dy.*dy.*eng); % no square ROOT; weighted
%     S(c) = sum(dy.*dy); % no square ROOT; not weighted
    dy = abs(y2-y1);
    S(c) = sqrt( sum(dy.*dy) );
end
score = - S;
W = ones(1,CN)./CN;

end

