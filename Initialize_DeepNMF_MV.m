function [Z, Hc, Hv] = Initialize_DeepNMF_MV(fea, layers, nviews)

    nOfLayers = length(layers);
    label = cell(nviews,1);

    Z = cell(nviews, nOfLayers);
    Hc = cell(1,nOfLayers);
    Hv = cell(nviews, nOfLayers);
    
    
    for layeri = 1:nOfLayers
        for t = 1:nviews
            rng('default');
            label{t} = litekmeans(fea{1,t}, layers(layeri));  % doc
            for i = 1:layers(layeri)
                Hv{t,layeri}(label{t} ==i, i) = 1;
            end

            Hv{t,layeri} = (Hv{t,layeri}+0.2)'; %H{view1,layer1} size 100 x nsamp
            mfea_viewt = size(fea{1,t},2);
            if layeri ==1
                Z{t,layeri} = ones(mfea_viewt, layers(layeri));
            else
                Z{t,layeri} = ones(layers(layeri-1), layers(layeri));
            end 
        end
        for i = 1:layers(layeri)
            Hc{1,layeri}(label{t} ==i, i) = 1;
        end
        Hc{1,layeri} = (Hc{1,layeri}+0.2)';
    end
    

    %{
    for t = 1:nviews
        for layeri = 1:nOfLayers%eg 2 layers 100 50
            rng('default');
            label{t} = litekmeans(fea{1,t}, layers(layeri));  % doc
            for i = 1:layers(layeri)
                if(t==2)
                    Hc{1,layeri}(label{t} ==i, i) = 1;
                end
                Hv{t,layeri}(label{t} ==i, i) = 1;
            end
            if(t==2)
                Hc{1,layeri} = (Hc{1,layeri}+0.2)';
            end
            Hv{t,layeri} = (Hv{t,layeri}+0.2)'; %H{view1,layer1} size 100 x nsamp
            mfea_viewt = size(fea{1,t},2);
            if layeri ==1
                Z{t,layeri} = ones(mfea_viewt, layers(layeri));
            else
                Z{t,layeri} = ones(layers(layeri-1), layers(layeri));
            end    
        end
    end
    %}

end

