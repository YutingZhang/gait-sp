function Info = zgResultArrStatistic( SArr  )


channelN = 5;
ccn = 0;
for c=1:channelN
    ccn = ccn + nchoosek( channelN, c );
end

InvalidArr = cellfun( @isempty, SArr );
ValidIndex = ~InvalidArr(:);

% --- Rank 1 Identification rates
Info.Rank1 = repmat( { NaN( size(SArr) ) }, ccn, 1 );
ValidSArr = cell2mat(SArr(ValidIndex));
IDENT = {ValidSArr.IDENT};
IDENT = cat( 2, IDENT{:} );
IDENT_SUM = cellfun( @sum, IDENT );

for cc = 1:ccn
    rc  = cat( 1, IDENT{cc,:} );
    r1c = rc(:,1)./IDENT_SUM(cc,:).';
    Info.Rank1{cc}(ValidIndex) = r1c;
end

% --- EER
Info.EER      = ArrangeEER( { ValidSArr.AUTH }, InvalidArr );
Info.EER_NORM = ArrangeEER( { ValidSArr.AUTH_NORM }, InvalidArr );


end


function EER = ArrangeEER( AUTH , InvalidArr )

ValidIndex = ~InvalidArr(:);
AUTH = cat( 2, AUTH{:} );
EERa = reshape([ AUTH.EER ],size(AUTH));

ccn = size(AUTH,1);
EER = repmat( { NaN( size(InvalidArr) ) }, ccn, 1 );
for cc = 1:ccn
    ec  = EERa(cc,:).';
    EER{cc}(ValidIndex) = ec;
end


end