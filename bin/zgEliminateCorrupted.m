function [ R W PersonNum ] = zgEliminateCorrupted( R, W, origPersonNum, Corrupted_Read, Corrupted_Test )
% [ R W ] = zgEliminateCorrupted( R, W, origPersonNum, Corrupted_Read, Corrupted_Test )

bothNum   = origPersonNum.both;
singleNum = origPersonNum.single;
totalNum  = bothNum+singleNum;

ReadIDX = true( totalNum, 1 );
ReadIDX(Corrupted_Read) = false;
TestIDX = true( totalNum, 1 );
TestIDX(Corrupted_Test) = false;

if ~isstruct(R)

    if ~all(TestIDX(ReadIDX))
        error( 'Required person is not in the test result' );
    end

    PIDX = ReadIDX(TestIDX);

    R = R(:,:,:,PIDX,PIDX,:);
    W = W(:,:,:,PIDX,PIDX,:);
    
else
    
    if ~all(ReadIDX == TestIDX)
        error( 'Two CorruptedPerson options are not consistent' );
    end
    
    
end

PersonNum.both   = sum(ReadIDX(1:bothNum));
PersonNum.single = sum(ReadIDX((bothNum+1):end));

end

