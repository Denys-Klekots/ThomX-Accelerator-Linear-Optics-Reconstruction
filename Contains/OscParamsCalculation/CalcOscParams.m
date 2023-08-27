function [Amplitudes, Phases] = CalcOscParams(A,s)
    % Calculate amplitudes and phases for one separate oscilation
    % A is the N_bpms * 2 matrix
    % s is the 2 * N_turns matrix
    % return vector columns on N_bpms length
    % phases restores only in range from  tp 2*pi and with starting
    % point of phases(0) = 0, as there is only phase difference has a sense 
    
    [~, N_turns] = size(s);
    [N_bpms, ~] = size(A);

    [cosMode, sinMode] = permutModes(s);
    
    Amplitudes = sqrt( A(:,1).^2 + A(:,2).^2 )*sqrt(2/N_turns);

    Phases = ones(N_bpms,1);
    for i=1:1:N_bpms
      tanPsiX = - A(i, sinMode)/A(i, cosMode);
      if abs(tanPsiX) < 1
        Phases(i) = atan(tanPsiX);
      else
        Phases(i) = acot(1/tanPsiX);
      end

      if A(i, cosMode) < 0
        Phases(i) = Phases(i) + pi;
      end

    end
    
 Phases = Phases - Phases(1);
 Phases = mod(Phases, 2*pi);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cosMode, sinMode] = permutModes(s)

  mode1s = s(1, :)/max(abs(s(1, :)));
  mode2s = s(2, :)/max(abs(s(2, :)));

  compValues = mode1s + 1i*mode2s;
  ang = angle(compValues);
  ang = ang(1:ceil(length(ang)/8));

  ang = ang(2:end) - ang(1:end-1);

  cos1Idx = (abs(ang) < pi) .* (ang > 0) + (abs(ang) > pi) .* (ang < 0);
  cos2Idx = (abs(ang) < pi) .* (ang < 0) + (abs(ang) > pi) .* (ang > 0);

  if sum(cos1Idx)/2 > sum(cos2Idx)*2
    cosMode = 1;
    sinMode = 2;
  elseif sum(cos1Idx)/2 < sum(cos2Idx)*2
    cosMode = 2;
    sinMode = 1;
  else
    ME = MException('CalcOcsParams:unableSeparateModes', ...
        'The number of signs of cos and sin mode is equal');
    throw(ME);    
  end

end