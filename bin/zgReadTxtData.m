function [BOTH_DATA SINGLE_DATA] = zgReadTxtData()

% INIT
recordCNT = 6;

src_path = '../data/txt2';

both_dir_file   = [ src_path '/record_both.dir' ];
single_dir_file = [ src_path '/record_single.dir' ];

both_list   = readLines( both_dir_file );
single_list = readLines( single_dir_file );

% BOTH
N = length(both_list);
BOTH_DATA = cell(recordCNT,N,2);
for k = 1:N
    for d = 1:2
        for r = 1:recordCNT
            record_name =  [ both_list{k} '-' int2str(d) '-' int2str(r) ];
            BOTH_DATA{r,k,d} = readRecord( src_path, record_name );
        end
    end
end

% SINGLE
N = length(single_list);
SINGLE_DATA = cell(recordCNT,N);
for k = 1:N
    for r = 1:recordCNT
        record_name =  [ single_list{k} '-' int2str(r) ];
        SINGLE_DATA{r,k} = readRecord( src_path, record_name );
    end
end

end

function T = readLines( file_path )

T = {};
k = 0;
fid = fopen( file_path, 'r' );
tline = fgetl(fid);
while ischar(tline)
    k = k+1;
    T{k}  = tline;
    tline = fgetl(fid);
end
fclose(fid);

end

function a = readRecord( src_path, record_name )

channelCNT = 5;

a.record_name = record_name;

spFN = [ src_path '/' a.record_name '.splitter' ];
if exist( spFN , 'file' )
    a.splitters = dlmread( spFN );
else
    fprintf(1,'Splitter missing: %s\n',spFN);
    a.splitters = [];
end

% a.data = dlmread( [ src_path '/' a.record_name '-1.txt' ] );
for c = 1:channelCNT
    a.data(:,:,c) = dlmread( [ src_path '/' a.record_name '-' int2str(c) '.txt' ] );
    cycFN = [ src_path '/' a.record_name '-' int2str(c) '.cycle' ];
    if exist( cycFN, 'file' )
        a.cycles{c} = dlmread( cycFN );
    else
        a.cycles{c} = [];
    end
end

if numel(a.splitters)<1
    a.splitters = [1 size(a.data,2)];
end
if numel(a.splitters)<2
    a.splitters(2) = size(a.data,2);
end

end