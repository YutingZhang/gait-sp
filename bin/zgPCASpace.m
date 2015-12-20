function [Dic Xm Eng] = zgPCASpace( X, energy_threshold )

if ~exist('energy_threshold','var')
    energy_threshold = 1;
end

Xm = mean(X,2);

Xc = X-repmat(Xm,1,size(X,2)); 

[U S V] = svd(Xc);

Eng = diag(S).^2;
EngS= sum(Eng);
if EngS>0
    Eng = Eng./EngS;
end
CEng = cumsum(Eng);
CEng(end) = 1;
idx  = find( CEng>=energy_threshold-eps, 1 );
Eng = Eng(1:idx);
Dic = U(:,1:idx);

end
