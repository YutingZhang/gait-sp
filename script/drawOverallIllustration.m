isOLD = 0;

clf

load ../results/spsrc/feat.mat
load gait.mat
load score_param.mat

availablePersonIdx = true( size( FEAT.BOTH, 2 ), 1 );
availablePersonIdx( score_param.CorruptedPerson ) = false;
FEAT.BOTH = FEAT.BOTH(:,availablePersonIdx,:);
DATA.BOTH = DATA.BOTH(:,availablePersonIdx,:);

c  = 5;
d1 = 1;
sr = 1;
er = 1;
sk = 1;
d2 = 3 - d1;

s  = 3;

load( sprintf( '../results/my-detailed/BaseInfo-c%d.mat', c ) );
load( sprintf( '../results/my-detailed/ident-c%d-d%d-sr%d-er%d-sk%03d.mat', ...
    c, d1, sr, er, sk ) );

[~, GLabel]= max( Rarr, [], 1 );

idxSelectedGSig = find( GLabel == sk, s );
idxSelectedGSig = idxSelectedGSig(end);

sGSigs = As{sr,d2}( :, Bs{sr,d2} == sr );
g = sGSigs( :, idxSelectedGSig );

chosenCycleID = FEAT.BOTH{sr,sk,d2}{c}.cycleID( idxSelectedGSig );
chosenLoc     = FEAT.BOTH{sr,sk,d2}{c}.loc( idxSelectedGSig );

stPos = DATA.BOTH{sr,sk,d2}.cycles{5}(chosenCycleID);
enPos = DATA.BOTH{sr,sk,d2}.cycles{5}(chosenCycleID+1);

D = DATA.BOTH{sr,sk,d2}.data(:,:,c);
D = sum( D.*D, 1 );

gLoc = stPos+(enPos-stPos)*chosenLoc;
gVal = interp1( 1:length(D), D , gLoc );

% figure

axes
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [25 50 50 170*4-50], 'XLim', [0 50], 'YLim', [0 170*4-50], ...
    'XTick', [], 'YTick', [], 'box', 'off', 'Color', 'none', ...
    'XColor', 'white', 'YColor', 'white' );

text( 10, 170*3+50, 'On the probe series', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 11, 'FontName', 'Times' );

text( 10, 170+270/2, 'In the chosen sub-dictionary', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 11, 'FontName', 'Times' );

text( 10, 50, 'For final code', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 11, 'FontName', 'Times' );

line([12,12],[170*3+50-70,170+270/2+110],[0 0], 'Color', 'Black');
line([12,12+3],[170+270/2+110,170+270/2+110+6],[0 0], 'Color', 'Black');
line([12,12-3],[170+270/2+110,170+270/2+110+6],[0 0], 'Color', 'Black');

line([12,12],[170+270/2-105,50+65],[0 0], 'Color', 'Black');
line([12,12+3],[50+65,50+65+6],[0 0], 'Color', 'Black');
line([12,12-3],[50+65,50+65+6],[0 0], 'Color', 'Black');


axes
plot(stPos:enPos, D(stPos:enPos) );
hold all
plot( gLoc, gVal, 'or', 'LineWidth', 1 );
plot( gLoc, gVal, '.g' );
hold off
set( gca, 'XLim', [stPos, enPos] );

set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100 50+170*3 250 100] );
xlabel( 'Time / ms' );
ylabel( 'Acceleration / g' );
title( 'a) A step cycle in the probe series' );

axes
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100 50+170*3 250 100], 'Color', 'none', ...
    'box', 'off', 'XLim', [0 1], 'YLim', [0 1], 'XTick', [], 'YTick', [] );
hold on
plot( 0.55, 0.85, 'or', 'LineWidth', 1 );
plot( 0.55, 0.85, '.g' );
text( 0.55, 0.85, '   Signature point', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle' );
hold off

axes
plot(g./norm(g));

set( gca, 'XLim',[ 1 21 ], 'XTick', 1:5:21, 'YLim', [0 0.5] );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100+300 50+170*3 150 100] );

xlabel( 'Dimension' );
title( 'b) SP descriptor' );

idxChosenClusters = find( chosenClusterIdxB(:,idxSelectedGSig) );

% figure
axes
plot( WordsIdent{er,d1}(:,idxChosenClusters) );

set( gca, 'XLim',[ 1 21 ], 'XTick', 1:5:21, 'YLim', [0 0.5] );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100+300 50+170*2 150 100] );

xlabel( 'Dimension' );

title( 'd) Sub-dictionary center' );

idxChosenDic = ismember( ClusterLabelIdent{er,d1}, idxChosenClusters ); 

chosenAe = Ae{er,d1}( :, idxChosenDic );
chosenBe = Be{er,d1}( idxChosenDic );

if isOLD
    idxChosenDic( length(Le{er,d1})+1:end ) = false; % patch, eliminate it later
end
chosenLe = Le{er,d1}( idxChosenDic );

axes
plot( chosenAe );

set( gca, 'XLim',[ 1 21 ], 'XTick', 1:5:21, 'YLim', [0 0.5] );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100+300 50+170*1 150 100] );

xlabel( 'Dimension' );

title( 'e) Sub-dictionary elements' );

% figure
axes
h = [];
h(1) = plot( chosenLe, -ones(size(chosenLe)), '.r' );
line( [0 1],[0 0], 'Color', 'black' );
hold on
[y x] = ksdensity( chosenLe );
h(2) = plot( x,y, '-b' );
hold off
set( gca, 'XLim',[ 0 1 ], 'YLim', [-2 ceil(max(y))], 'YTick', 0:2:ceil(max(y)) );
set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100 50+170*1 250 270] );
set( h(1), 'DisplayName', 'Point locations' );
set( h(2), 'DisplayName', 'Point density' );
hl = legend( h, 'Location', 'Best' );
set( hl, 'EdgeColor', 'white' );
xlabel( 'In-cycle relative location (gait phase)' );
ylabel( 'Density' );

title( 'c) Locational distribution for sub-dictionary' );

% final coding

% if isOLD
%     ga = omp( chosenAe, g, [], 8 );
%     ga = full(ga);
%     ga_parted  = data2cellByIndex( ga , chosenBe, 1, size(Rarr,1) );
%     gan1 = cell( @(x) sum(abs(gan1)), ga_parted );
% else
    gan1 = GaNorm1( :, idxSelectedGSig );
% end

axes
h = [];
h(2) = stem( gan1([1:sk-1,sk+1:end]), 'b', 'MarkerSize', 4 );
hold on
h(1) = stem( gan1(sk), 'r' );
hold off
set( gca, 'XLim', [0 176], 'XTick', [ 0:10:length(gan1) ] );

set( gca, 'Unit', 'pixel' );
set( gca, 'Position', [100 50+0 450 100] );

xlabel( 'Subject Id.' );
ylabel( 'Coefficient l1-norm' );
title( 'f) Coefficient distribution with respect to subjects' );

set( h(1), 'DisplayName', [ 'Subject ' int2str(sk) ] );
set( h(2), 'DisplayName', 'Other subjects' );
hl = legend( h, 'Location', 'Best' );
set( hl, 'EdgeColor', 'white' );


