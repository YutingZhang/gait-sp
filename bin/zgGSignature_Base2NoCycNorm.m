function D = zgGSignature_Base2NoCycNorm( BASE_D, RECORD )

CN = length( BASE_D );
D=cell(1,CN);

for c=1:CN
    A = BASE_D{c};
    AIDX = and(A.loc>=RECORD.splitters(1),A.loc<=RECORD.splitters(2));
    A.descriptor = A.descriptor(AIDX,:);
    A.loc   = A.loc(AIDX);
    A.scale = A.scale(AIDX);
    D{c} = A;
end

end
