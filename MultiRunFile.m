function MultiRunFile(datasetname, weightMode, normalizationType)
    
    % Writing all the outputs to a text file including errors (if occured),
    % any print text on console and etc.
    %fileName = strcat(extractBefore(datasetname,'.mat'),'_',normalizationType,'_',weightMode,'.txt');
    %diary (fileName);
    
    addpath(genpath('.'));
    warning off;
    
    fprintf('Dataset name: %s / Normalization Type: %s / Weight Mode: %s \n\n',datasetname, normalizationType, weightMode);
    load(append('data/',datasetname));
    
    temp = cell(1,numel(fea));
    for i = 1:nviews
        temp{1,i} = NormalizeFea(transpose(fea{i,1}),0);
    end 
    CorrectFea = temp;
    CorrectGnd = gnd;

    %%% Start From Here %%%%

    nClass = length(unique(CorrectGnd)); %document cluster number
    nviews = numel(CorrectFea); %number of views

    originalDataMatrix = cell(1,nviews);

    for i = 1:nviews
        %Taking transpose since row wise concatination is required
        originalDataMatrix{1,i} = transpose(CorrectFea{1,i});
    end

    rng('default');

    %First view matrix or concatenate all views to be clustered
    %firstView = originalDataMatrix{1,1};
    %concat = horzcat(originalDataMatrix{:});
    
    %{
    %Clustering in the original space
    disp('Clustering in the original space, concatenating of all views - Kmeans.');
    [best.CA, best.F, best.P, best.R, best.nmi, best.AR] = performance_kmeans(concat, nClass, CorrectGnd);
    disp(['    NMI and std:       ',num2str(best.nmi(1)), ' , ', num2str(best.nmi(2))]);
    disp(['    Accuracy and std:  ',num2str(best.CA(1)), ' , ', num2str(best.CA(2))]);
    disp(['    F-score and std:   ',num2str(best.F(1)), ' , ', num2str(best.F(2))]);
    %}

    layers = [100 50];
    k = 30;
    alpha = 100;
    alphaForL = 1;
    beta = 100;
    delta = 0.001;
    maxIteration = 200;
    expRawCount = 1;
    currentCount = 1;
    totalMetrics = 9;

    % This can have two values #1 denotes non-negative constraints and #2 denotes mix signs. 
    % However, we are only using #1 at the moment 
    algorithm_no = 1;
    % Set this flag to 0 if parameter training is not required
    isParaTrainReq = 0;

    if isParaTrainReq ==1 
        
        layersArray = [100 50 0 0; 150 50 0 0; 150 100 50 0; 100 75 50 0; 200 150 100 50];
        kArray = [20 30 50 70 100 120 130 140 150 180];
        alphaArray = [0.1 0.5 1 10 20 50 80 100];
        betaArray = [0.001 0.01 0.1 0.5 1 5 10 20 50 100];
        deltaArray = [0.001 0.01 0.1 0.5 1 5 10 20 50 100];
        
        expRawCount = getTotalRawCount(layersArray, kArray, alphaArray, betaArray, deltaArray);
        currentCount = 0;

        HcKmeansResults = cell(expRawCount,totalMetrics);
        HcSpectralResults = cell(expRawCount,totalMetrics); 
        HdKmeansResults = cell(expRawCount,totalMetrics); 
        HdSpectralResults = cell(expRawCount,totalMetrics); 
        AllKmeansResults = cell(expRawCount,totalMetrics); 
        AllSpectralResults = cell(expRawCount,totalMetrics);
        
        for layersArrayIndex = 1 : length(layersArray)

            layer1 = layersArray(layersArrayIndex,1);
            layer2 = layersArray(layersArrayIndex,2);
            layer3 = layersArray(layersArrayIndex,3);
            layer4 = layersArray(layersArrayIndex,4);

            if(layer3 ==0 && layer4 ==0)
                layers = [layer1 layer2];
            elseif(layer4 ==0)
                layers = [layer1 layer2 layer3];
            else
                layers = [layer1 layer2 layer3 layer4];
            end

            for kArrayIndex = 1 : length(kArray)
                k = kArray(kArrayIndex);
    
                for alphaArrayIndex = 1 : length(alphaArray)
                    alpha = alphaArray(alphaArrayIndex);
    
                    for betaArrayIndex = 1 : length(betaArray)
                        beta = betaArray(betaArrayIndex);
    
                        for deltaArrayIndex = 1 : length(deltaArray)
                            delta = deltaArray(deltaArrayIndex);
                            
                            currentCount = currentCount + 1;
                            [HcKmeansResults, HcSpectralResults, HdKmeansResults, HdSpectralResults, AllKmeansResults, AllSpectralResults] = DoClustering(layers, k, alpha, beta, delta, maxIteration, originalDataMatrix, nviews, alphaForL, algorithm_no, CorrectGnd, weightMode, nClass, datasetname, currentCount, HcKmeansResults, HcSpectralResults, HdKmeansResults, HdSpectralResults, AllKmeansResults, AllSpectralResults, normalizationType);
                        end
                    end
                end
            end
        end

    else
        HcKmeansResults = cell(expRawCount,totalMetrics);
        HcSpectralResults = cell(expRawCount,totalMetrics); 
        HdKmeansResults = cell(expRawCount,totalMetrics); 
        HdSpectralResults = cell(expRawCount,totalMetrics); 
        AllKmeansResults = cell(expRawCount,totalMetrics); 
        AllSpectralResults = cell(expRawCount,totalMetrics);
        
        DoClustering(layers, k, alpha, beta, delta, maxIteration, originalDataMatrix, nviews, alphaForL, algorithm_no, CorrectGnd, weightMode, nClass, datasetname, currentCount, HcKmeansResults, HcSpectralResults, HdKmeansResults, HdSpectralResults, AllKmeansResults, AllSpectralResults, normalizationType);
    end

    %diary off

    return
end

function [HcKmeansResults, HcSpectralResults, HdKmeansResults, HdSpectralResults, AllKmeansResults, AllSpectralResults] = DoClustering(layers, k, alpha, beta, delta, maxIteration, originalDataMatrix, nviews, alphaForL, algorithm_no, CorrectGnd, weightMode, nClass, datasetname, currentCount, HcKmeansResults, HcSpectralResults, HdKmeansResults, HdSpectralResults, AllKmeansResults, AllSpectralResults, normalizationType)
    
    fprintf('\n Clustering output when iter = %f / layers = %s / k = %f / alpha = %f / beta = %f / delta = %f \n\n', maxIteration, layers, k, alpha, beta, delta);
    
    [ Z, Hc, Hv, maxItrError ] = DCCNMF_Function(originalDataMatrix, k, alpha, beta, delta, nviews, layers, maxIteration, alphaForL, algorithm_no, weightMode, normalizationType);

    if(~isnan(maxItrError) && ~isinf(maxItrError))
        nOfLayers = length(layers);

        %FinalHv = horzcat(transpose(Hv{1,nOfLayers}), transpose(Hv{2,nOfLayers}), transpose(Hv{3,nOfLayers}));
        
        noOfRowsinHv = size(transpose(Hv{1,nOfLayers}),1);
        noOfColumnsinHv = size(transpose(Hv{1,nOfLayers}),2);

        FinalHv = zeros(noOfRowsinHv,noOfColumnsinHv,nviews) ; 
        for i = 1:nviews
            FinalHv(:,:,i) = transpose(Hv{i,nOfLayers}) ; 
        end

        FinalHv = reshape(FinalHv,noOfRowsinHv,noOfColumnsinHv*nviews);

        %Hd Spectral
        disp('Clustering using concatination for all Hd - Spectral');
        [best.CA, best.F, best.P, best.R, best.nmi, best.AR] = performance_SpectralClustering(transpose(FinalHv), CorrectGnd);
        disp(['    NMI and std:       ',num2str(best.nmi(1)), ' , ', num2str(best.nmi(2))]);
        disp(['    Accuracy and std:  ',num2str(best.CA(1)), ' , ', num2str(best.CA(2))]);
        disp(['    F-score and std:   ',num2str(best.F(1)), ' , ', num2str(best.F(2))]);
        disp(['    ARI and std:       ',num2str(best.AR(1)), ' , ', num2str(best.AR(2))]);

        HvNMI = best.nmi(1);
        HvAcc = best.CA(1);
        HvFscore = best.F(1);
        HvARI = best.AR(1);

        HdSpectralResults{currentCount,1} = HvNMI;
        HdSpectralResults{currentCount,2} = HvAcc;
        HdSpectralResults{currentCount,3} = HvFscore;
        HdSpectralResults{currentCount,4} = HvARI;
        
        HdSpectralResults{currentCount,5} = layers;
        HdSpectralResults{currentCount,6} = k;
        HdSpectralResults{currentCount,7} = alpha;
        HdSpectralResults{currentCount,8} = beta;
        HdSpectralResults{currentCount,9} = delta;

        HvFileName = strcat('HdResults_Spectral_',weightMode,'_',normalizationType,'_',datasetname);
        save(HvFileName,'HdSpectralResults');

    end
end

function totalCount = getTotalRawCount(layersArray, kArray, alphaArray, betaArray, deltaArray)
    
    totalCount = 0;
    for layersArrayIndex = 1 : length(layersArray)
        for kArrayIndex = 1 : length(kArray)
            for alphaArrayIndex = 1 : length(alphaArray)
                for betaArrayIndex = 1 : length(betaArray)
                    for deltaArrayIndex = 1 : length(deltaArray)
                        totalCount = totalCount + 1;
                    end
                end
            end
        end
    end
end