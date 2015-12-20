function D = zgGSignature_Base2WithCycNorm( BASE_D, RECORD )

%% Call non-normalized version
D = zgGSignature_Base( RECORD );

%% Processing
CN = length( BASE_D );
cycles = RECORD.cycles{5};  % defautly use channel 5
cycleLEN = diff(cycles);



[ cycST cycEN ] = zgCycleEnds( cycles , RECORD.splitters );


for c=1:CN
    A   = D{c};
    loc = D{c}.loc;
    IDX  = quantiz(loc,cycles);
    AIDX = and(IDX>=cycST,IDX<cycEN);
    IDX  = IDX(AIDX);
    A.descriptor = A.descriptor(AIDX,:);
    A.scale = A.scale(AIDX);
    loc = loc(AIDX);
    A.loc = (loc - cycles(IDX))./cycleLEN(IDX);
    A.loc(A.loc == 1) = 0;
    A.baseCycleID = cycST;
    A.cycleID = IDX - (cycST-1);
    D{c} = A;
end

end

