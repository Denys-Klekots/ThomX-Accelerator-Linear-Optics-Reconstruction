function DrawLinearParameters(AM, imgFoldername)

    mkdir(imgFoldername)
    

    ring = GetLattice(AM);
    inex_BPM = find(atgetcells(ring, 'Class', 'Monitor'))';
    
    [refLinParam, refT, refC] = atlinopt(ring, 0, inex_BPM);
    
    betaXRef = cat(1, refLinParam.beta);
    betaXRef = betaXRef(:,1);
    
    betaYRef = cat(1, refLinParam.beta);
    betaYRef = betaYRef(:,2);

    refPsiX = cat(1, refLinParam.mu);
    refPsiX = refPsiX(:,1);
    refPsiX = refPsiX-refPsiX(1);

    refPsiY = cat(1, refLinParam.mu);
    refPsiY = refPsiY(:,2);
    refPsiY = refPsiY-refPsiY(1);

    
    dispRef = cat(1, refLinParam.Dispersion);
    dispRef = dispRef(1:4:end);



    
    [x_readings, y_readings] = ExtractBPMReadings(AM);

    N_samples = length(x_readings);

    [N_bpms, ~] = size(x_readings{1});

    DispNorm = [];
    DispMean = [];

    BetaXNorm = [];
    BetaXMean = [];

    BetaYNorm = [];
    BetaYMean = [];

    PhiD = [];
    PhiX = [];
    PhiY = [];

    SynchFreq = [];
    BetaXFreq = [];
    BetaYFreq = [];

    for smpN = 1:1:N_samples

        X = x_readings{smpN};
        X = X(:,20:420);
        X = X - mean(X,2);

        Y = y_readings{smpN};
        Y = Y(:,20:420);
        Y = Y - mean(Y,2);

        [L, A, s] = ICA([X;Y], N_bpms);

        try
%             [AmpS, PhS, freqS] = CalcSampleParamsFreqBounds(A, ... %boundaries for callibration
%                                                        s, 0.00, 0.25, L);
%             [AmpX, PhX, freqX] = CalcSampleParamsFreqBounds(A, ...
%                                                        s, 0.32, 0.50, L);
%             [AmpY, PhY, freqY] = CalcSampleParamsFreqBounds(A, ...
%                                                        s, 0.25, 0.31, L);


%             [AmpS, PhS, freqS] = CalcSampleParamsFreqBounds(A, ... %boundaries for files with y oscilations
%                                                          s, 0.00, 0.25, L);
%             [AmpX, PhX, freqX] = CalcSampleParamsFreqBounds(A, ...
%                                                          s, 0.35, 0.50, L);
%             [AmpY, PhY, freqY] = CalcSampleParamsFreqBounds(A, ...
%                                                          s, 0.25, 0.35, L);


            [AmpS, PhS, freqS] = CalcSampleParamsFreqBounds(A, ... %boundaries for files with only x oscilations
                                                         s, 0.00, 0.25, L);
            [AmpX, PhX, freqX] = CalcSampleParamsFreqBounds(A, ...
                                                         s, 0.26, 0.50, L);
            AmpY = zeros(2*N_bpms,1);% y plane is unneeded in this case
            PhY  = zeros(2*N_bpms,1);
            freqY = 0;


            % correcting values, considering the features of BPM readings
            freqX = freqX/2;
            AmpX = AmpX(1:N_bpms) ./ cos(pi*freqX);
            PhX = PhX(1:N_bpms);     
            AmpX = AmpX/1000000000; %translate amplitude in meter units


            freqS = freqS/2;
            AmpS = AmpS(1:N_bpms) ./ cos(pi*freqS);
            PhS = PhS(1:N_bpms);     
            AmpS = AmpS/1000000000; %translate amplitude in meter units            


            %correcting do tune that nu_y > 0.5 and feature of BPM 
            freqY = (1+freqY)/2;
            AmpY = AmpY(N_bpms+1:end) ./ cos(pi*freqY);
            % PhY need additional mod 2*pi as the PhY(0) = 0,
            % not PhY(N_bpms+1) = 0
            PhY = mod(PhY(N_bpms+1:end) - PhY(N_bpms+1), 2*pi); 
            AmpY = AmpY/1000000000; %translate amplitude in meter units            


        catch ME
            if (strcmp(ME.identifier,'CalcOcsParams:unableSeparateModes'))
                disp("Unable to distribute modes, skipping sample")
                continue
            else
                rethrow(ME);
            end

        end     

        DispMean = [DispMean, mean(AmpS, 1)];
        DispNorm = [DispNorm, AmpS/DispMean(end)];        


        BetaXMean = [BetaXMean, mean(AmpX.^2, 1)];
        BetaXNorm = [BetaXNorm, (AmpX.^2)/BetaXMean(end)];

        BetaYMean = [BetaYMean, mean(AmpY.^2, 1)];
        BetaYNorm = [BetaYNorm, (AmpY.^2)/BetaYMean(end)];        

        PhiD = [PhiD, PhS .* zeros(N_bpms,1)];
        PhiX = [PhiX, PhX];
        PhiY = [PhiY, PhY];

        SynchFreq = [SynchFreq, freqS];
        BetaXFreq = [BetaXFreq, freqX];
        BetaYFreq = [BetaYFreq, freqY];

    end

    [~,goodSamplesN] = size(SynchFreq);
    if goodSamplesN <3
     return
    end


    DrawFigures(BetaXMean, BetaXNorm, betaXRef, PhiX, refPsiX,BetaXFreq,...
                "Horisontal Plane", "\beta_x [m]", "\psi_x [rad]",...
                "\nu_x [turn^{-1}]","\epsilon_x [m rad]", imgFoldername);
        
    DrawFigures(DispMean, DispNorm, dispRef, PhiD, zeros(size(PhiD)), ...
                SynchFreq, "Longitudinal motion", " D [m]", ...
                "\psi_s [rad]", "\nu_s [turn^{-1}]", ...
                " (\Delta p/p_0)_{max}", imgFoldername);

%     DrawFigures(BetaYMean, BetaYNorm, betaYRef, PhiY, refPsiY,BetaYFreq,...
%                 "Vertical Plane", "\beta_y [m]", "\psi_y [rad]",...
%                 "\nu_y [turn^{-1}]","\epsilon_y [m rad]", imgFoldername);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DrawFigures(AMean, ANorm, ARef, Ph, PhRef, Freq,...
              Title, Alabel, PhLabel, FreqLabel, ScaleLabel, imgFoldername)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Amplitude values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    AScale = AMean/mean(ARef);

    [N_bpms, N_samples] = size(ANorm);

    ANorm = ANorm*mean(ARef);    
    ANormMean = mean(ANorm,2)';
    
    histA = reshape(ANorm', N_bpms*N_samples,1);

    xHist = ceil( (1:1:N_samples*N_bpms)/N_samples) - 0.25;
    xHist = xHist';

    xages = 0.5:0.5:N_bpms+0.5;
    yages = linspace(0, max(histA), 25);
    
    
    fig=figure;
    subplot(2,1,1);
    
    histogram2(xHist, histA, "XBinEdges", xages, ...
               "YBinEdges", yages, 'FaceColor','flat', 'DisplayStyle',...
               'tile','ShowEmptyBins','off', "EdgeColor", "none");
    colormap(flipud(hot));
    colorbar;
    
    % displaying amplitude value restored from BPM and reference
    hold on;
    plot([1:1:N_bpms]+0.25, ANormMean, "b-o");
    hold on;
    plot([1:1:N_bpms]+0.25, ARef, "g-o");
    title(Title);
    xlabel("BPM index");
    ylabel(Alabel);

    % displaying the amplitude value, restored from BPMs
    % and the reference in separate plot
    subplot(2,1,2);
    plot(ANormMean', "b-o");
    hold on;
    plot(ARef, "g-o");
    title(Title);
    xlabel("BPM index");
    ylabel(Alabel);
    legend("Restored", "Simulated", "Location", 'eastoutside');  

    Alabel = regexprep(Alabel, "[\\\/]", "");
    figName = strcat(imgFoldername, "/",Title, "_", Alabel, ".png");
    saveas(fig, figName)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Phase values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PhMean = exp(1i*Ph);
    PhMean = mean(PhMean, 2);
    PhMean = angle(PhMean);

    for i=2:1:N_bpms
        while PhMean(i) < PhMean(i-1)
            PhMean(i) = PhMean(i) + 2*pi;
        end
    end

    Ph = Ph + floor( (PhMean - Ph)/2/pi )*2*pi;
    Ph( abs(PhMean - Ph) > abs(PhMean - Ph - 2*pi) ) = ...
                   Ph( abs(PhMean - Ph) > abs(PhMean - Ph - 2*pi) ) + 2*pi;


    histPh = reshape(Ph', N_bpms*N_samples,1);

    yages = linspace(0, max(histPh), 25);

    fig = figure;
    subplot(2,1,1);
    
    histogram2(xHist, histPh, "XBinEdges", xages, ...
               "YBinEdges", yages, 'FaceColor','flat', 'DisplayStyle',...
               'tile','ShowEmptyBins','off', "EdgeColor", "none");
    colormap(flipud(hot));
    colorbar;    

    % displaying phase value restored from BPM and reference    
    hold on;
    plot([1:1:N_bpms]+0.25, PhMean, "b-o");
    hold on;
    plot([1:1:N_bpms]+0.25, PhRef, "g-o");
    title(Title);
    xlabel("BPM index");
    ylabel(PhLabel);


    % displaying the phase vale, restored from BPMs
    % and the reference in separate plot
    subplot(2,1,2);
    plot(PhMean', "b-o");
    hold on;
    plot(PhRef, "g-o");
    title(Title);
    xlabel("BPM index");
    ylabel(PhLabel);
    legend("Restored", "Simulated", "Location", 'eastoutside');  


    PhLabel = regexprep(PhLabel, "[\\\/]", "");
    figName = strcat(imgFoldername, "/",Title, "_", PhLabel, ".png");
    saveas(fig, figName)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Scale and Frequency values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    fig = figure;
    histogram(AScale, linspace(min(AScale)-1E-10, max(AScale)+1E-10, 10));
    title(Title)
    xlabel(ScaleLabel);

    ScaleLabel = regexprep(ScaleLabel, "[\\\/]", "");    
    figName = strcat(imgFoldername, "/",Title, "_", ScaleLabel, ".png");
    saveas(fig, figName)


    fig = figure;
    histogram(Freq, linspace(min(Freq)-0.01, max(Freq)+0.01, 10));
    title(Title)
    xlabel(FreqLabel);

    FreqLabel = regexprep(FreqLabel, "[\\\/]", "");      
    figName = strcat(imgFoldername, "/",Title, "_", FreqLabel, ".png");
    saveas(fig, figName)

    close all
end
