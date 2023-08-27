function [Amp, Ph, freq, modeNums] = ...
          CalcSampleParamsFreqBounds(A, s, lBound, rBound, L)
% calculate parameters of the oscilations in one data sample
% A is spatial matrix from ICA with dimensions N_bpms * P
% s is the templorary matrix crom ICA with dimentions P*N_turns
% lBound and rBound is the bound in which the oscilation will be searched,
%        the function will work with the two first modes, which has the
%        largest peack in range between boundaries. 
% L is the square of the eigen value of the ICA analysis
%
% Amp is the N_bpms length column vector of the amplitudes of oscilations 
% Ph similar to the Amp but contain phases
% freq contain frequency of the oscilation
% mode nums is the 2 element array, containing the number os used modes

    [N_bpms,~] = size(A);
    Amp = zeros(N_bpms, 1);
    Ph = zeros(N_bpms, 1);

    [modeNums, freq] =  DistributeModes(s.*sqrt(L), lBound, rBound);

    [Amp(:), Ph(:)] = CalcOscParams(A(:,modeNums),s(modeNums,:));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [modeNums, freq] =  DistributeModes(Ls, lBound, rBound)
% Ls is the s matrix from ICA, but multiplied on eigen values of the ICA
% analysis
%
% modeNums is the 2 element array containing mode numbers
% freq is the freauensy of oscilation    
    
    [N_modes, N_turns] = size(Ls);

    lBoundIndex = floor(lBound*N_turns);
    rBoundIndex = ceil(rBound*N_turns);

    peakFreqs = ones(1, N_modes);
    peakValues = ones(1, N_modes);

    for modN=1:1:N_modes
        spectr = abs( fft(Ls(modN,:)) );
        spectr = spectr(1:1:floor(N_turns/2));

        spectX = (1:1:floor(N_turns/2))/N_turns;

        [pks,loc] = findpeaks(spectr, spectX,'SortStr',...
                              'descend', 'NPeaks', 1);

        peakFreqs(modN) = loc(1);

        if loc(1) > lBound && loc(1) < rBound
            peakValues(modN) = pks(1);
        else
            peakValues(modN) = 0;
        end

    end

    

    [maxPeakValues(1), modeNums(1)] = max(peakValues);
    peakValues(modeNums(1)) = 0;
    [maxPeakValues(2), modeNums(2)] = max(peakValues);

    % weighet mean value
    freq = (peakFreqs(modeNums(1))*maxPeakValues(1) + ...
            peakFreqs(modeNums(2))*maxPeakValues(2))/...
            (maxPeakValues(1) + maxPeakValues(2));


end


