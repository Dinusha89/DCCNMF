% Normalization is performed in column wise since samples are represented
% using column
function Matrix = NormalizeHv(Matrix, nOfLayers, nviews, typeOfNormalization)      
    if(typeOfNormalization == "None")
        return;
    else    
        for i_layer = 1:nOfLayers
            for v_ind = 1:nviews
                if(typeOfNormalization == "L1")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 1, 'norm',1);
                elseif(typeOfNormalization == "L2")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 1, 'norm',2);
                elseif(typeOfNormalization == "MinMax")
                    Matrix{v_ind,i_layer} = normalize(Matrix{v_ind,i_layer}, 1, 'range');
                elseif(typeOfNormalization == "Std")
                    Matrix{v_ind,i_layer} = zscore(Matrix{v_ind,i_layer}, 0, 1);
                end
            end
        end
    end
end