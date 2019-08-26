function Results = Predict_pairs(encoded_training_alignment, encoded_focus_alignment, LengthA, table_count_species)
%makes pairing predictions on encoded_focus_alignment using Mirrortree

[N, alignment_width] = size(encoded_focus_alignment);
L=alignment_width-2; % last 2 columns contain species index and initial sequence index

%initialize the Results array, used for saving data
Results=zeros(N-2,5);
%col 1: species
%col 2: HK index in initial alignment
%col 3: RR index in initial alignment
%col 4: score of pairing
%col 5: gap 

%total pair counter
totcount = 0; 


%loop over species
for i=1:size(table_count_species,1)
    
    test_seqs = encoded_focus_alignment(table_count_species(i,2):table_count_species(i,3),:);
    NSeqs = table_count_species(i,3)-table_count_species(i,2)+1;
    species_id=table_count_species(i,1);
    
    distA=pdist2(encoded_training_alignment(:,1:LengthA),test_seqs(:,1:LengthA),'hamm');
    %a column contains distances of one sequence A of test_seqs to those of the training set
    distB=pdist2(encoded_training_alignment(:,LengthA+1:L),test_seqs(:,LengthA+1:L),'hamm');

    %now compute the correlations of all the distance vectors
    HKRR_energy = -corr(distA,distB);
    %and compute the predicted assignment
    if NSeqs==1
        assignment=1;
        HKRR_energy_b=HKRR_energy-min(min(HKRR_energy)); %ensure that all elements are >=0
    elseif isequal(min(HKRR_energy(:)),max(HKRR_energy(:)))
        assignment=randperm(NSeqs); %avoids spurious positive results
        HKRR_energy_b=HKRR_energy-min(min(HKRR_energy)); %ensure that all elements are >=0
    else %use the Hungarian algorithm
        HKRR_energy_b=HKRR_energy-min(min(HKRR_energy)); %ensure that all elements are >=0
        [assignment, score] = assignmentoptimal(HKRR_energy_b);
        %deal with identical rows, i.e. effectively identical HKs
        uEn = unique(HKRR_energy_b, 'rows');
        if size(uEn,1) < size(HKRR_energy_b,1)
            assignment=randomize_equal_rows(assignment,HKRR_energy_b,uEn);
        end
        %deal with identical cols, i.e. effectively identical RRs
        uEn = unique(HKRR_energy_b', 'rows'); %transpose to deal with columns
        if size(uEn,1) < size(HKRR_energy_b,2) %uEn is transposed
            assignment=randomize_equal_cols(assignment,HKRR_energy_b,uEn);
        end
    end
    
    bigval=1e3*abs(max(HKRR_energy_b(:)));
    
    for j=1:NSeqs
        totcount=totcount+1;
        Results(totcount,1)=species_id;
        Results(totcount,2)=test_seqs(j,L+2);  %initial index of HK sequence (HK: line)
        Results(totcount,3)=test_seqs(assignment(j),L+2); %initial index of RR sequence (RR: col)
        Results(totcount,4)=HKRR_energy(j,assignment(j)); %absolute energy of the pairing
        if NSeqs==1
            Results(totcount,5)=abs(HKRR_energy); %no real gap... consider that absolute energy is gap
        elseif isequal(min(HKRR_energy(:)),max(HKRR_energy(:)))
            %no gap for this assignment
            Results(totcount,5)=0;
        else
            %calculate gap for this assignment
            HKRR_energy_mod=HKRR_energy_b;
            HKRR_energy_mod(j,assignment(j))=bigval;
            [~, score_mod] = assignmentoptimal(HKRR_energy_mod);
            Results(totcount,5)=score_mod-score;
        end
    end
    
end

end
