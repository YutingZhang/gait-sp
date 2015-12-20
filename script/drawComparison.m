needSave2pdf = 1;

load results

close all

CC = [ 3 25 31 ];

lColor = lines(8);

% CPlotStyleBase = { ...
%     { 'Color', 'red' } ; % Cycle matching
%     { 'Color', [1 0.5 0] } ; % Average cycle
%     { 'Color', 'blue' } ; % Eigen step
%     { 'Color', 'green' } ; % Self-DTW4Gallery
%     { 'Color', 'green' } ; % SP voting
%     { 'Color', 'green' } ; % msSP voting
%     { 'Color', 'black' }; };

CPlotStyleBase = { ...
    { 'Color', lColor(1,:) } ; % Cycle matching
    { 'Color', lColor(2,:) } ; % Average cycle
    { 'Color', lColor(3,:) } ; % Eigen step
    { 'Color', lColor(4,:) } ; % Self-DTW4Gallery
    { 'Color', lColor(5,:) } ; % SP voting
    { 'Color', lColor(7,:) } ; % msSP voting
    { 'Color', 'black' }; };

CPlotStyle = { ...
    { '-.' } ; % Cycle matching
    { '--' } ; % Average cycle
    { '-.' } ; % Eigen step
    { '--' } ; % Self-DTW4Gallery
    { '-.' } ; % SP voting
    { '--' } ; % msSP voting
    { '-' }; };

CMarkerStyle = { ...
    { 'Marker', '+', 'MarkerSize', 5 } ; % Cycle matching 
    { 'Marker', 'x', 'MarkerSize', 6 } ; % Average cycle
    { 'Marker', 's', 'MarkerSize', 5 } ; % Eigen step
    { 'Marker', 'd', 'MarkerSize', 6 } ; % Self-DTW4Gallery
    { 'Marker', 'v', 'MarkerSize', 5 } ; % SP voting
    { 'Marker', '^', 'MarkerSize', 5 } ; % msSP voting
    { 'Marker', 'o', 'MarkerFace', 'green', 'MarkerSize', 5 } };
    

% CMC

A = [ 1 1 1 0 1 1 1 ];
lengthCMC = 10;

for r = 1:length(CC)
    figure( 'name', [ 'CMC: ' int2str(CC(r))] );
    set( gca, 'YLim', [0 100], 'YTick', 0:10:100, ...
        'XLim', [1-0.5 lengthCMC+0.5], 'XTick', 1:1:lengthCMC, ...
        'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
    set( gca, 'Unit', 'pixel' );
    set( gca, 'Position', [ 50 50 350 350 ] );
    hold on;
    h = []; kk = 0;
    for k = 1:length(CArr)
        if A(k)
            CMC_k = cumsum( CArr{k}.IDENT{CC(r)} )./sum(CArr{k}.IDENT{CC(r)});
            CMC_k = CMC_k * 100;
            kk = kk+1;
            h(kk) = plot( 1:lengthCMC, CMC_k(1:lengthCMC), ...
                CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
            set( h(kk), 'DisplayName', CName{k} );
        end
    end
    xlabel( 'Rank' );
    ylabel( 'Identification rate / %' );
    legend(h, 'Location', 'SouthEast' );
    hold off;
    if needSave2pdf
        savepdf( ['CMC-' int2str(CC(r))] );
    end
end

% Location number - IR

BN = [{1:5},{6:15},{16:25},{26:30},{31}];



figure( 'name', [ 'IR: BN' ] );
set( gca, 'YLim', [0 100], 'YTick', 0:10:100, ...
    'XLim', [0.5 length(BN)+0.5], 'XTick', 1:1:length(BN), ...
    'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [ 50 50 350 350 ] );
hold on;
h = []; kk = 0;
for k = 1:length(CArr)
    if A(k)
        r1BN = cell2mat( CInfo{k}.Rank1 );
        averageR1BN = cellfun( @(x) mean(r1BN(x)), BN);
        kk = kk+1;
        h(kk) = plot( 1:1:length(BN), averageR1BN*100, ...
            CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
        set( h(kk), 'DisplayName', CName{k} );
    end
end
xlabel( 'Number of body locations' );
ylabel( '(Mean) identification rate / %' );
legend(h, 'Location', 'SouthEast' );
hold off;
if needSave2pdf
    savepdf( ['IR-BN'] );
end

% ROC

A = [ 1 1 1 0 0 0 1 ];

hasUnnorm = [ 1 1 1 0 0 0 0 ];
hasUnnorm(:) = 0;

for r = 1:length(CC);
    figure( 'name', [ 'ROC: ' int2str(CC(r))] );
    set( gca, ...
        'YLim', [0 100], 'YTick', 0:10:100, ...
        'XLim', [0 100], 'XTick', 0:10:100, ...
        'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
    set( gca, 'Unit', 'pixel' );
    set( gca, 'Position', [ 50 50 350 350 ] );
    hold on;
    
    h = []; kk = 0;
    for k = 1:length(CArr)
        if A(k)
            TPN  = CArr{k}.AUTH_NORM(CC(r)).ACCEPT_TOTAL - CArr{k}.AUTH_NORM(CC(r)).FALSE_REJECT;
            FPR  = CArr{k}.AUTH_NORM(CC(r)).FALSE_ACCEPT / CArr{k}.AUTH_NORM(CC(r)).REJECT_TOTAL;
            TPR  = TPN / CArr{k}.AUTH_NORM(CC(r)).ACCEPT_TOTAL;    % recall (sensitivity)
            [ uFPR uidx ] = unique(FPR);
            uTPR = TPR(uidx);
            markerFPR = 0:0.05:1;
            markerFPR = markerFPR(2:end-1);
            markerTPR = interp1( uFPR, uTPR, markerFPR, 'linear' );
            plot( FPR*100, TPR*100, CPlotStyle{k}{:}, CPlotStyleBase{k}{:} );
            plot( markerFPR*100, markerTPR*100, ...
                CMarkerStyle{k}{:}, CPlotStyleBase{k}{:}, 'LineStyle', 'none' );

            kk = kk+1;
            h(kk) = plot( -1, -1, ...
                CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
            set( h(kk), 'DisplayName', CName{k} );

            if hasUnnorm(k)
                FPR = CArr{k}.AUTH(CC(r)).FALSE_ACCEPT / CArr{k}.AUTH(CC(r)).REJECT_TOTAL;
                TPR = 1 - CArr{k}.AUTH(CC(r)).FALSE_REJECT / CArr{k}.AUTH(CC(r)).ACCEPT_TOTAL;
                kk    = kk+1;
                h(kk) = plot( FPR*100, TPR*100, CPlotStyle{k}{:}, CPlotStyleBase{k}{:} );
                set( h(kk), 'DisplayName', [ CName{k} ' - raw' ] );
            end
        end
        
    end
    xlabel( 'False positive rate / %' );
    ylabel( 'True positive rate / %' );
    legend(h, 'Location', 'SouthWest' );
    
    hold off;
    
    if needSave2pdf
        savepdf( ['ROC-' int2str(CC(r))] );
    end
    
    % ---------------------------------------------------
    
    figure( 'name', [ 'PR: ' int2str(CC(r))] );
    set( gca, ...
        'YLim', [0 100], 'YTick', 0:10:100, ...
        'XLim', [0 100], 'XTick', 0:10:100, ...
        'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
    set( gca, 'Unit', 'pixel' );
    set( gca, 'Position', [ 50 50 350 350 ] );
    hold on;
    
    h = []; kk = 0;
    for k = 1:length(CArr)
        if A(k)
            TPN  = CArr{k}.AUTH_NORM(CC(r)).ACCEPT_TOTAL - CArr{k}.AUTH_NORM(CC(r)).FALSE_REJECT;
            PREC = TPN ./ ( CArr{k}.AUTH_NORM(CC(r)).FALSE_ACCEPT + TPN ); % precision
            RECA = TPN ./   CArr{k}.AUTH_NORM(CC(r)).ACCEPT_TOTAL;         % recall (sensitivity)
            [ uRECA uidx ] = unique(RECA);
            uPREC = PREC(uidx);
            markerRECA = 0:0.05:1;
            markerRECA = markerRECA(2:end-1);
            markerPREC = interp1( uRECA, uPREC, markerRECA, 'linear' );
            plot( RECA*100, PREC*100, CPlotStyle{k}{:}, CPlotStyleBase{k}{:} );
            plot( markerRECA*100, markerPREC*100, ...
                CMarkerStyle{k}{:}, CPlotStyleBase{k}{:}, 'LineStyle', 'none' );

            kk = kk+1;
            h(kk) = plot( -1, -1, ...
                CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
            set( h(kk), 'DisplayName', CName{k} );

            if hasUnnorm(k)
                TPN   = CArr{k}.AUTH(CC(r)).ACCEPT_TOTAL - CArr{k}.AUTH(CC(r)).FALSE_REJECT;
                PREC  = TPN ./ ( CArr{k}.AUTH(CC(r)).FALSE_ACCEPT + TPN ); % precision
                RECA  = TPN ./   CArr{k}.AUTH(CC(r)).ACCEPT_TOTAL;    % recall (sensitivity)
                kk    = kk+1;
                h(kk) = plot( RECA*100, PREC*100, CPlotStyle{k}{:}, CPlotStyleBase{k}{:} );
                set( h(kk), 'DisplayName', [ CName{k} ' - raw' ] );
            end
        end
        
    end
    xlabel( 'Recall / %' );
    ylabel( 'Precision / %' );
    legend(h, 'Location', 'SouthEast' );
    
    hold off;
    
    if needSave2pdf
        savepdf( ['PR-' int2str(CC(r))] );
    end
    
end

% Location number - EER

BN = [{1:5},{6:15},{16:25},{26:30},{31}];

figure( 'name', [ 'EER: BN' ] );
set( gca, 'YLim', [0 100], 'YTick', 0:10:100, ...
    'XLim', [0.5 length(BN)+0.5], 'XTick', 1:1:length(BN), ...
    'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [ 50 50 350 350 ] );
hold on;
h = []; kk = 0;
for k = 1:length(CArr)
    if A(k)
        r1BN = cell2mat( CInfo{k}.EER_NORM );
        averageR1BN = cellfun( @(x) mean(r1BN(x)), BN);
        kk = kk+1;
        h(kk) = plot( 1:1:length(BN), averageR1BN*100, ...
            CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
        set( h(kk), 'DisplayName', CName{k} );
    end
end
xlabel( 'Number of body locations' );
ylabel( '(Mean) equal error rate / %' );
legend(h, 'Location', 'NorthEast' );
hold off;
if needSave2pdf
    savepdf( ['EER-BN'] );
end


% CYCLE

load cycle_results

A = [ 1 1 1 0 1 1 1 ];

ccTitle = {'GP', 'P', 'G'};
ccType  = {'(gallery & probe)', 'probe', 'gallery'};

for r = 1:length(CC)
    for j = 1:size(CArrCyc,1)
        
        figure( 'name', [ 'Controlled Cycles - ' ccTitle{j} ' : ' int2str(CC(r))] );
        set( gca, ...
            'YLim', [0 100], 'YTick', 0:10:100, ...
            'XLim', [1-0.5 10+0.5], 'XTick', 1:1:10, ...
            'box', 'on', 'XGrid', 'off', 'YGrid', 'on' );
        set( gca, 'Unit', 'pixel' );
        set( gca, 'Position', [ 50 50 350 350 ] );
        hold on;

        h = []; kk = 0;
        for k = 1:size(CArr,2)
            if A(k)
                IR = squeeze( CInfoCyc{j,k}.Rank1{CC(r)} );
                kk = kk + 1;
                h(kk) = plot( 1:length(IR), IR*100, ...
                    CPlotStyle{k}{:}, CMarkerStyle{k}{:}, CPlotStyleBase{k}{:} );
                set( h(kk), 'DisplayName', CName{k} );
            end
        end
        
        xlabel( [ 'Cycle number per ' ccType{j} ' series for' ] );
        ylabel( 'Identification rate / %' );
        legend(h, 'Location', 'SouthEast' );
        hold off;
        
        if needSave2pdf
            savepdf( ['CC-' ccTitle{j} '-' int2str(CC(r))] );
        end
    end
end


