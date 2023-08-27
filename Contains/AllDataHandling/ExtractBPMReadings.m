function [x_readings, y_readings] = ExtractBPMReadings(AM)

    N_samples = length(AM.Data.RingBPM);
    N_bpms = length(AM.Data.RingBPM(1).dev_replies);
    N_tunrs = length(AM.Data.RingBPM(1).dev_replies(1 ...
                     ).attr_values(3).value);

    x_readings = cell(1,N_samples);
    y_readings = cell(1,N_samples);

    for readingNum = 1:1:N_samples
       X = ones(N_bpms, N_tunrs);
       Y = ones(N_bpms, N_tunrs);

       for bpmNum = 1:1:N_bpms
          X(bpmNum,:) = AM.Data.RingBPM(readingNum).dev_replies(bpmNum...
                        ).attr_values(3).value;

          Y(bpmNum,:) = AM.Data.RingBPM(readingNum).dev_replies(bpmNum...
                        ).attr_values(4).value;

       end

        permutations = 2:1:N_bpms;
        permutations = [permutations, 1];
        X = X(permutations,:);
        Y = Y(permutations,:);

       x_readings(readingNum) = {X};
       y_readings(readingNum) = {Y};
    end
end

