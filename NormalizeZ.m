% Normalization is performed in row wise since samples are represented
% using rows
function Matrix = NormalizeZ(Matrix, nOfLayers, nviews, typeOfNormalization)      
    if(typeOfNormalization == "None")
        return;
    else    
        for i_layer = 1:nOfLayers
            for v_ind = 1:nviews
                if(typeOfNormalization == "L1")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 2, 'norm',1);
                elseif(typeOfNormalization == "L2")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 2, 'norm',2);
                elseif(typeOfNormalization == "MinMax")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 2, 'range');
                elseif(typeOfNormalization == "Std")
                    Matrix{v_ind,i_layer} = zscore(Matrix{v_ind,i_layer}, 0, 2);
                end
            end
        end
    end
end