function [descriptor loc scale] = zgSIFT1D(seq, dspHalfSize, Smax, sigma0, sigmaN, layer_per_oct, threshold)
% loc = SIFTFind1D(seq, Smax, [,sigma0[, sigmaN [,threshold]]])
% loc{Smax+1}(keypoint_num)
% Smin = 0

%% Handling args
if nargin<3
    Smax = 2;
end
if nargin<4
    sigma0 = 1.6;
end
if nargin<5
    sigmaN = 0.5;
end
if nargin<6
    layer_per_oct = 2;  % scale layer per oct
end
if nargin<7
    threshold = 0.01/layer_per_oct/sigma0;
end

if(Smax<0)
    error('Smax must not be less than 0');
end
if(sigma0<=0)
    error('sigma0 must be greater than 0');
end
if(sigmaN<0)
    error('sigmaN must not be less than 0');
end
if(layer_per_oct < 1)
    error('layer_per_oct must be greater or equal to 1');
end

%% Initialize
% 

scaleNum = layer_per_oct * Smax;

seqLen = length(seq);

gau = zeros( scaleNum+3, seqLen );
% DoG = zeros( scaleNum+2, seqLen );
extrema = zeros( scaleNum+2, seqLen  );

% 1 & 2 order differentials
Dd1x = zeros(scaleNum+2,seqLen);
Dd2x = zeros(scaleNum+2,seqLen);
Dd1s = zeros(scaleNum+2,seqLen);
Dd2s = zeros(scaleNum+2,seqLen);

%% Generate Guassians

alpha = 2^(1/layer_per_oct);


for m = 0:scaleNum+2
    sigma = sqrt((sigma0 * alpha^(-2+m))^2 - (sigmaN)^2);
    gau(m+1,:) = normfilter( seq, sigma );
end

DoG = diff( gau, 1, 1 );

Dd1x(:,2:end-1) = ( DoG(:,3:end) - DoG(:,1:end-2) )/2;
Dd2x(:,2:end-1) = ( DoG(:,3:end) + DoG(:,1:end-2) - 2*DoG(:,2:end-1) );
Dd1s(2:end-1,:) = ( DoG(3:end,:) - DoG(1:end-2,:))/2;
Dd2s(2:end-1,:) = ( DoG(3:end,:) + DoG(1:end-2,:) - 2*DoG(2:end-1,:) );

%% find extrema & calculate differentials


Neigh = cat( 3 ...
    , DoG(1:end-2,2:end-1) ...  % 1
    , DoG(1:end-2,1:end-2) ...  % 2
    , DoG(2:end-1,1:end-2) ...  % 3
    , DoG(3:end  ,1:end-2) ...  % 4
    , DoG(3:end  ,2:end-1) ...  % 5
    , DoG(3:end  ,3:end  ) ...  % 6
    , DoG(2:end-1,3:end  ) ...  % 7
    , DoG(1:end-2,3:end  ) ...  % 8
    );

Neigh = bsxfun( @minus, Neigh, DoG(2:end-1,2:end-1) );
% --- ^^^^ efficiency can be improved by 2 times

extrema(2:end-1,2:end-1) = ( abs( sum(sign(Neigh),3) ) >= 8 );

%% refine location & scale

MaxIterNum = 3;

[ M N ] = find( extrema );

%% offset

I = sub2ind( size(DoG), M, N );
I = reshape( I, numel(I), 1 );
workingIdx = true( size(I) );
X_off = zeros( size(N) );
S_off = zeros( size(M) );

% for iter=1:MaxIterNum;
% 
%     if iter>1
%         N(workingIdx) = N(workingIdx) + round( X_off( workingIdx ) );
%         M(workingIdx) = M(workingIdx) + round( S_off( workingIdx ) );
%         Mcmin = and( M<1,workingIdx );  Mcmax = and( M>size(DoG,1),workingIdx );
%         Ncmin = and( N<1,workingIdx );  Ncmax = and( N>size(DoG,2),workingIdx );
%         M( Mcmin ) = 1; M( Mcmax ) = size(DoG,1);
%         N( Ncmin ) = 1; N( Ncmax ) = size(DoG,2);
%         workingIdx( Mcmin | Mcmax | Ncmin | Ncmax ) = false;
%         if ~any(workingIdx)
%             break;
%         end
%         I(workingIdx) = sub2ind( size(DoG), M(workingIdx), N(workingIdx) );
%     end
%     
%     J = I(workingIdx);
%     
%     X_off(workingIdx) = -Dd1x(J)./(Dd2x(J)+1e-5);
%     S_off(workingIdx) = -Dd1s(J)./(Dd2s(J)+1e-5);
% 
%     % rollback unstable estimation
%     X_off(X_off>=1) = 0;
%     S_off(S_off>=1) = 0;
% 
%     % update workingIdx
%     workingIdx( and( X_off<=0.5, S_off<=0.5 ) ) = false;
% 
%     % stop criterion
%     if ~any(workingIdx)
%         break;
%     end
% 
%     
% end

% rule out unstable ones

T1 = [Dd1x(I),Dd1s(I)].*[X_off,S_off];
validKP = ( abs( DoG(I)+0.5*(T1(:,1)+T1(:,2)) ) > threshold );

S = M(validKP) + S_off(validKP);
X = N(validKP) + X_off(validKP);


%% Compute descriptors


scale = (alpha.^(-2+S)) * sigma0;
validKP = and( (-dspHalfSize*scale+X)>=1, (dspHalfSize*scale+X)<=seqLen );

S     = S(validKP);
X     = X(validKP);

scale       = scale(validKP);
loc         = X;
descriptor  = zeros( length(S), dspHalfSize*2+1 );

Sp = S; Sp(Sp<1) = 1;
sampleLayer = repmat( reshape(Sp,numel(Sp),1), 1, size(descriptor,2) );
sampleIDX   = bsxfun( @plus, reshape(X,numel(X),1), ...
    bsxfun( @times, -dspHalfSize:dspHalfSize, reshape(scale,numel(scale),1) ) );

descriptor(:) = interp2( gau, sampleIDX(:), sampleLayer(:) );


end

%% normfilter
function fld = normfilter(sq,sigma)

    if sigma <= 0
        fld = sq;
    elseif sigma > 0
        bound = ceil(sigma*3);
        fg = normpdf(-bound:bound, 0, sigma);
        fg = fg/sum(fg);
        fld = convfilter(fg, sq);
    else
        error('sigma must not be less than 0');
    end

end

%% filter

function y = convfilter( kernel, x )
    y = conv( x, kernel(end:-1:1), 'same' );
end