function [CA, F, P, Recall, nmi, AR] = performance_SpectralClustering( H, truth)

    nClass = length(unique(truth));
    replic = 10; % Number of replications for Spectral Clustering

    if iscell(H)
        H = H{:};
    end

    for i = 1:replic

        options = [];
        options.k = 5;
        options.WeightMode = 'HeatKernel';
        options.t = 1;
        A = constructW(H', options);
        gt = truth; clusNum = nClass;

        idx = SpectralClustering(A,clusNum);

        %dgunlk
        %{
        [A nmii(i) avgent] = compute_nmi(gt,C);
        [F,P,R] = compute_f(gt,C);
        [AR,RI,MI,HI]=RandIndex(gt,C);
        C = bestMap(gt,C);
        ACC(i) = length(find(gt == C))/length(gt);
        %}
        
        CAi(i) = 1-compute_CE(idx, truth); % clustering accuracy
        [Fi(i),Pi(i),Ri(i)] = compute_f(truth,idx); % F1, precision, recall
        %dgunlk
        %[A, nmii(i), avgent] = compute_nmi(truth,idx);
        nmii(i) = compute_nmi(truth,idx);
        ARi(i) = rand_index(truth,idx);

    end

    CA(1) = mean(CAi); CA(2) = std(CAi);
    F(1) = mean(Fi); F(2) = std(Fi);
    P(1) = mean(Pi); P(2) = std(Pi);
    Recall(1) = mean(Ri); Recall(2) = std(Ri);
    nmi(1) = mean(nmii); nmi(2) = std(nmii);
    AR(1) = mean(ARi); AR(2) = std(ARi);

end

