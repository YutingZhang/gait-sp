function [R W paramPersonNum] = zgCrossValidationS( FEAT , score_fun , param, options )
% R = zgCrossValidationS( FEAT , score_fun , param )
if ischar(FEAT)
    load( FEAT );
end

paramPersonNum.both  = size(FEAT.BOTH,2);
paramPersonNum.single = size(FEAT.SINGLE,2);

person_total_both = size(FEAT.BOTH,2);
record_total = size(FEAT.BOTH,1);
person_total_single = size(FEAT.SINGLE,2);

% channel, sample record number, exemplar record number,
% sample persons, exemplar persons, days
R = zeros( [ 5, record_total , record_total , ...
    (person_total_both+person_total_single) , (person_total_both+person_total_single), 2 ] );
tS = size(R);
W  = zeros(tS);

SINGLE_D = FEAT.SINGLE;
for d1 = 1:2
    d2 = 3-d1;
    S_BOTH_D = FEAT.BOTH( :,:,d2 );
    for er=1:record_total
        E_BOTH   = FEAT.BOTH(er,:,d1);
        parfor ek = 1:person_total_both
            ts = tS;
            R0 = zeros([ts(1:2),1,ts(4)]);
            W0 = zeros( size(R0) );
            exemplar = E_BOTH{1,ek};
            S_BOTH   = S_BOTH_D;
            SINGLE   = SINGLE_D;
            for sr=1:record_total
                for sk=1:person_total_both
                    sample = S_BOTH{sr,sk};
                    [ R0( :, sr, 1, sk ), W0( :, sr, 1, sk ) ] = feval( score_fun, sample, exemplar, param );
                end
                for sk=1:person_total_single
                    sample = SINGLE{sr,sk};
                    [ R0( :, sr, 1, sk+person_total_both ), W0( :, sr, 1, sk+person_total_both ) ] = feval( score_fun, sample, exemplar, param );
                end
            end
            R(:,:,er,:, ek, d1) = R0;
            W(:,:,er,:, ek, d1) = W0;
        end
        clear('E_BOTH');
        parfor ek = 1:person_total_single
            ts = tS;
            R0 = zeros([ts(1:2),1,ts(4)]);
            W0 = zeros( size(R0) );
            exemplar = SINGLE_D{er,ek};
            S_BOTH   = S_BOTH_D;
            for sr=1:record_total
                for sk=1:person_total_both
                    sample = S_BOTH{sr,sk};
                    [ R0( :, sr, 1, sk ), W0( :, sr, 1, sk ) ] = feval( score_fun, sample, exemplar, param );
                end
            end
            R(:, :,er,:, ek+person_total_both, d1) = R0;
            W(:, :,er,:, ek+person_total_both, d1) = W0;
        end
        fprintf(1,'.');
    end
end
fprintf(1,'\n');

end
