function [Amp, Ph, freq] = CalcSampleParams(A, s, freqGuess)
    % calculate parameters of the oscilations in one data sample
    % A is spatial matrix from ICA with dimensions N_bpms * P
    % s is the templorary matrix crom ICA with dimentions P*N_turns
    % freqGues is the N_modes length row vector and contain initial guess
    %          of the frequency of each mode we interesting in
    %          note that function works only with first 2*N_modes columns
    %          of the A matrix and first 2*N_modes rows of s matrix
    %
    % Amp is the N_bpms * N_modes matrix, each column of it contain the
    %          the amplitudes of oscilations coresponding to freqGuess
    % Ph similar to the Amp but contain phases
    % freq coresponds to freqGues but contain more precise values

    N_osc = length(freqGuess);

    [N_bpms,~] = size(A);
    Amp = zeros(N_bpms, N_osc);
    Ph = zeros(N_bpms, N_osc);

    [modeNums, freq] =  DistributeModes(s, freqGuess);

    for oscN = 1:1:N_osc
        [Amp(:,oscN), Ph(:,oscN)] = ...
            CalcOscParams(A(:,modeNums(:,oscN)),s(modeNums(:,oscN),:));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [modeNums, freq] =  DistributeModes(s, freqGuess)
    % modeNums is the 2*length(freqGuess) mantrix, in columns contain 
    %       index of cos and sin modes of oscilations coresponding to
    %       freqGuess
    % freq is the row coresponding to dreqGuess but with precise values


    indx2Mode = ones(1,length(freqGuess));

    modeNums = zeros(2, length(freqGuess));

    freq = zeros(1, length(freqGuess));

    for modN = 1:1:2*length(freqGuess)

        spectr = abs( fft(s(modN,:)) );
        spectr = spectr(1:1:ceil(length(spectr)/2));

        spectX = ((1:1:length(spectr) )-1)/length(spectr)/2;       

        [~,loc] = findpeaks(spectr, spectX,'SortStr',...
                              'descend', 'NPeaks', 1);

        freqDev = abs( freqGuess - loc(1) );
        [~, freqIdx] = min(freqDev);
        
        modeNums(indx2Mode(freqIdx), freqIdx) = modN;
        indx2Mode(freqIdx) = indx2Mode(freqIdx) + 1;
        freq(freqIdx) = 0.5*loc(1) + freq(freqIdx);

        if indx2Mode(freqIdx) > 3
            ME = MException('CalcSampleParams:unableDistributeModes', ...
                'The freqGuess %i got third mode', freqGuess(modN));
            throw(ME);
        end

    end
end
