function drawLocationalRelevance( locationalHist )

clf

locationalHist = [ locationalHist(1:end-2), locationalHist(end-1)+locationalHist(end) ];
locationalHist = locationalHist./sum(locationalHist);

stepLength = 0.5/length( locationalHist );

bar( stepLength/2:stepLength:0.5, locationalHist );

set( gca, 'XLim', [-stepLength/2, 0.5+stepLength/2] , 'XTick', 0:stepLength:0.5 );

xlabel( ['Average distance to sub-dictionray / in-cycle relative location'] );
ylabel( 'Frequency' );

set( gca, 'Unit', 'Pixel' );
set( gca, 'Position', [100 50 350 200] );

end

