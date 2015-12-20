
if ~exist('DATA','var')
    load ../script/gait.mat
end


idSubj = num2cell([1,1,1]);
idCh   = 5;
idCycle= 1:3;

N = length(idCycle);

d = DATA.BOTH{idSubj{:}};
a = sqrt( sum( d.data(:,:,idCh).*d.data(:,:,idCh), 1 ));
c = d.cycles{5};

clf
for k=1:N
    subplot(N,1,k);
    cid = idCycle(k);
    plot( 0:1/(c(cid+1)-c(cid)):1, a( c(cid):c(cid+1) ) );
    
end
