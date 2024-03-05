% Normalization is performed in column wise since samples are represented
% using column
function Matrix = NormalizeHc(Matrix, nOfLayers, typeOfNormalization)
    if(typeOfNormalization == "None")
        return;
    else
        for i_layer = 1:nOfLayers
            if(typeOfNormalization == "L1")
                Matrix{1,i_layer} = normalize(Matrix{1,i_layer},1,'norm',1);
            elseif(typeOfNormalization == "L2")
                Matrix{1,i_layer} = normalize(Matrix{1,i_layer},1,'norm',2);
            elseif(typeOfNormalization == "MinMax")
                Matrix{1,i_layer} = normalize(Matrix{1,i_layer},1,'range');
            elseif(typeOfNormalization == "Std")
                Matrix{1,i_layer} = zscore(Matrix{1,i_layer}, 0, 1);
            end
        end
    end
end