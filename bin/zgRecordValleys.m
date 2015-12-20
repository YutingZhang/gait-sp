function [VAL MOD] = RecordValleys(ORIG, Channels, line_offset)
%[VAL MOD] = RecordValleys(ORIG [, Channels])

if ~exist('line_offset','var')
    line_offset=0;
end

%% Calculate Norms
DIS = zeros([size(ORIG,2) size(ORIG,3)]);
for k = 1:size(ORIG,3)
    for n = 1:size(ORIG,2)
        DIS(n,k) = norm(ORIG(:,n,k));
    end
end

peaks = zeros([100 size(DIS,2)]);

if nargin < 2
    Channels = 1:size(DIS,2);
end

%% find peaks for each channel
for k = Channels
    %% find raw minima ranges
    USEFUL_IDX = floor(size(DIS,1)*0.2):floor(size(DIS,1)*0.8);
    
    if 0
        pl = min(DIS(USEFUL_IDX,k));
        dc = mean(DIS(:,k));
        lt = pl + (dc - pl) * (0.45+line_offset);
    else
        sortedDIS1 = sort(DIS(USEFUL_IDX,k));
        lt = sortedDIS1(ceil(length(sortedDIS1)*(0.1+line_offset)));
    end
    
    if 0    % Only for debug
        fprintf( 1, '%f%%\n', sum(DIS(USEFUL_IDX,k)<lt)/length(DIS(USEFUL_IDX,k))*100 );
    end
    
    latent = DIS(:,k) < lt;
    %% find minima ranges
    ind = 1;
    ranges=zeros([100 2]);
    n = 1;
    while ind<=size(DIS,1)
        while (ind<=size(DIS,1) && latent(ind)==0)
            ind = ind+1;
        end
        rt = ind;
        while (ind<=size(DIS,1) && latent(ind)==1)
            ind = ind+1;
        end
        re = ind - 1;
        if (rt+3>re || re >= size(DIS,1)-1) % epsilon == 3
            continue;
        end
        ranges(n,:) = [rt re];
        n = n+1;
    end
    ranges=ranges(1:(n-1),:);
    %% find peaks with elimination of close peaks
    old_peak=1;
    n=1;
    for ind = 1:size(ranges,1)
        [cmin imin] = min(DIS(ranges(ind,1):ranges(ind,2),k));
        rindex = imin+ranges(ind,1)-1;
        if ind>1
%             if max(DIS(old_peak:rindex,k)) < lt
            if rindex-old_peak < 15
                if DIS(old_peak,k)>DIS(rindex,k)
                    n = n-1;
                else
                    continue;
                end
            end
        end
        peaks(n,k) = rindex;
        old_peak = peaks(n,k);
        n=n+1;
    end
    peak_num = n-1;
    %% add lost peaks
    tpeaks = peaks(1:peak_num,k);
    interval = tpeaks(2:peak_num) - tpeaks(1:peak_num-1);
    st_len = median(interval);
    if (st_len<50)
        st_len = 100;
    elseif (st_len<70)
        st_len = 70;
    elseif (st_len>140)
        st_len = 140;
    end
    intlarge = interval > (st_len * 1.7);
    ind = 1;
    for n = 1:peak_num-1
        peaks(ind, k) = tpeaks(n);
        ind=ind+1;
        if intlarge(n)
            split_num = round(interval(n) / st_len);
            split_interval = ceil(interval(n) / split_num);
            stp = tpeaks(n);
            stp = stp + split_interval;
            while stp<tpeaks(n+1)
                lowi = stp - ceil(st_len * 0.25);
                [cmin imin] = min(DIS(lowi:stp+ceil(st_len * 0.25),k));
                peaks(ind, k) = lowi+imin-1;
                ind = ind+1;
                stp = stp + split_interval;
            end
        end
    end
    peaks(ind, k) = tpeaks(peak_num);
    peak_num = ind;
    
    %% eliminate bad peaks
    tpeaks = peaks(1:peak_num,k);
    interval = tpeaks(2:peak_num) - tpeaks(1:peak_num-1);
    st_len = median(interval);
    n = 1;
    while n < peak_num
        ie = tpeaks(n)+st_len - st_len*0.3;
        m = n+1;
        while m<peak_num && tpeaks(m) < ie
            m = m+1;
        end
        if abs(tpeaks(m) - (tpeaks(n)+st_len)) < st_len*0.3
            break;
        end
        n=n+1;
    end
    ind = 1;
    while n <= peak_num
        peaks(ind,k)=tpeaks(n);
        ind = ind+1;
        old_peak = tpeaks(n);
        n = n+1;
        while n<=peak_num && tpeaks(n)-old_peak<=st_len*0.7
            n = n+1;
        end
        if(n==peak_num&&tpeaks(n)-old_peak<=st_len*0.7)
            break;
        end
    end
    for n=ind:peak_num
        peaks(n,k) = 0;
    end
end

%% set return value
VAL = cell([size(peaks,2) 1]);
for k = 1:size(peaks,2)
    VAL{k} = peaks(1:find(peaks(:,k),1,'last'),k);
end
MOD = DIS;
