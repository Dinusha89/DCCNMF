function [Z, Hc, Hv, maxItrError] = DCCNMF_Function(X, k, alpha, beta, delta, nviews, layers, maxIteration, alphaForL, algorithm_no, weightMode, NormalizationType)

    Av = cell(nviews,1);
    Avd = cell(nviews,1);
    Lvd = cell(nviews,1);
    
    options = [];
    options.k = k;
    options.WeightMode = weightMode;
    options.NeighborMode = 'KNN';
    options.NormW = 1;
    
    nOfLayers = length(layers);
    %ErrorToPlot = zeros(maxIteration);
    
    for v_ind = 1:nviews
        %Taking transpose to present data in expected manne
        if algorithm_no == 1
            R = abs(transpose(X{1,v_ind}));
        elseif algorithm_no == 2
            R = transpose(X{1,v_ind});
        end
        
        %R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
        % Equation 1
        Av{v_ind,1} = constructW(R', options);
        Av{v_ind,1} = full(Av{v_ind,1});
        
    end
    for v_ind = 1:nviews
        Ac(:,:,v_ind) = Av{v_ind,1};
    end
    
    % Equation 2
    Ac = min(Ac,[],3);
    
    for v_ind = 1:nviews
        % Equation 3
        Avd{v_ind,1} = Av{v_ind,1} - Ac;
        % Equation 5
        Lvd{v_ind,1} = constructL(Avd{v_ind,1}, alphaForL, options);
    end
    
    % Equation 4
    Lc = constructL(Ac, alphaForL, options);
    
    LcPlus = (abs(Lc) + Lc)/2;
    LcMinus = (abs(Lc) - Lc)/2;
    
    LvdPlus = cell(nviews,1);
    LvdMinus = cell(nviews,1);
    
    [Z, Hc, Hv] = Initialize_DeepNMF_MV(X, layers, nviews);
    
    Hc = NormalizeHc(Hc, nOfLayers, NormalizationType);
    Hv = NormalizeHv(Hv, nOfLayers, nviews, NormalizationType);
    Z = NormalizeZ(Z, nOfLayers, nviews, NormalizationType);
    
    nIteration = 0;
    %CurrentDiffError = 1;
    
    while(maxIteration > nIteration)
        nIteration = nIteration + 1;
        
        for i_layer = 1:nOfLayers

            for v_ind = 1:nviews
                
                if algorithm_no == 1
                    R = abs(transpose(X{1,v_ind}));
                elseif algorithm_no == 2
                    R = transpose(X{1,v_ind});
                end
                %R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));

                LvdPlus{v_ind,1} = (abs(Lvd{v_ind,1}) + Lvd{v_ind,1})/2;
                LvdMinus{v_ind,1} = (abs(Lvd{v_ind,1}) - Lvd{v_ind,1})/2;
                
                % update Zvi / Equation 18 and 19
                Z = updateZ(R, Z, Hc, Hv, v_ind, i_layer, algorithm_no);
                Z = NormalizeZ(Z, nOfLayers, nviews, NormalizationType);
                
                % update Hvd / Equation 12
                Hv = updateHv(R, Z, Hc, Hv, beta, v_ind, i_layer, LvdPlus, LvdMinus);
                Hv = NormalizeHv(Hv, nOfLayers, nviews, NormalizationType);

            end
            % update Hc / Equation 9
            Hc = updateHc(X, Z, Hc, Hv, nviews, alpha, delta, i_layer, LcPlus, LcMinus, algorithm_no);
            Hc = NormalizeHc(Hc, nOfLayers, NormalizationType);
            
        end        
        
        % calculate error using the cost function / Equation 6
        CurrentError = cost_function(X, Z, Hc, Hv, Lc, Lvd, nviews, nOfLayers, alpha, beta, delta);
        
        %{
        if(nIteration>1)
            ErrorToPlot(nIteration) = CurrentError;
        end
        fprintf('Iteration = %d / Error = %f \n',nIteration, CurrentError);
        %}
        
        
        %Printing error only in last iteration for parameter tuning
        if(maxIteration == nIteration)
            fprintf('\n Iteration = %d / Error = %f \n\n',nIteration, CurrentError);
            maxItrError = CurrentError;
        end
        
    end
    %plot(ErrorToPlot);
    
end

function Z = updateZ(R, Z, Hc, Hv, v_ind, i_layer, algorithm_no)
    
    if algorithm_no == 1
        if i_layer == 1
            DX = R;
            DD = Z{v_ind, i_layer};                    
        else
            DX = CalculateDXZ(Z, R, v_ind, i_layer);
            DD = CalculateDDZ(Z, v_ind, i_layer);
        end

        DXPlus = (abs(DX) + DX)/2;
        DXMinus = (abs(DX) - DX)/2;

        DDPlus = (abs(DD) + DD)/2;
        DDMinus = (abs(DD) - DD)/2;

        % Equation 16                
        Zup = DXPlus*(transpose(Hc{1,i_layer}) + transpose(Hv{v_ind, i_layer})) + DDMinus*(Hc{1,i_layer}*transpose(Hc{1,i_layer}) + Hc{1,i_layer}*transpose(Hv{v_ind, i_layer}) + Hv{v_ind, i_layer}*transpose(Hc{1,i_layer}) + Hv{v_ind, i_layer}*transpose(Hv{v_ind, i_layer}));

        % Equation 17
        Zun = DXMinus*(transpose(Hc{1,i_layer}) + transpose(Hv{v_ind, i_layer})) + DDPlus*(Hc{1,i_layer}*transpose(Hc{1,i_layer}) + Hc{1,i_layer}*transpose(Hv{v_ind, i_layer}) + Hv{v_ind, i_layer}*transpose(Hc{1,i_layer}) + Hv{v_ind, i_layer}*transpose(Hv{v_ind, i_layer}));

        % Equation 18        
        Z{v_ind, i_layer} = Z{v_ind, i_layer}.*power((Zup./Zun),(0.5));
    
    elseif algorithm_no == 2
        if i_layer == 1
            Zup = R * (transpose(Hc{1,i_layer}) + transpose(Hv{v_ind, i_layer}));
            Zun = Z{v_ind, i_layer}*(Hc{1,i_layer}*transpose(Hc{1,i_layer}) + Hc{1,i_layer}*transpose(Hv{v_ind, i_layer}) + Hv{v_ind, i_layer}*transpose(Hc{1,i_layer}) + Hv{v_ind, i_layer}*transpose(Hv{v_ind, i_layer}) );
        else
            ThetaTransposeWithMinus = CalculateThetaTransposeWithMinus(Z, v_ind, i_layer);
            Theta = CalculateTheta(Z, v_ind, i_layer);
            
            Zup = ThetaTransposeWithMinus * R * ( transpose(Hc{1,i_layer}) + transpose(Hv{v_ind, i_layer}) );
            Zun = ThetaTransposeWithMinus * Theta * ( Hc{1,i_layer}*transpose(Hc{1,i_layer}) + Hc{1,i_layer}*transpose(Hv{v_ind, i_layer}) + Hv{v_ind, i_layer}*transpose(Hc{1,i_layer}) + Hv{v_ind, i_layer}*transpose(Hv{v_ind, i_layer}));
        end

        % Equation 19
        Z{v_ind, i_layer} = Z{v_ind, i_layer}.*(Zup./max(Zun, 1e-10));
    end

    Z{v_ind, i_layer}(isnan(Z{v_ind, i_layer}))=0;

end

function Hc = updateHc(X, Z, Hc, Hv, nviews, alpha, delta, i_layer, LcPlus, LcMinus, algorithm_no)
        
    % verify how multi-view handles in this scenario
    DX = CalculateDXHc(Z, X, nviews, i_layer, algorithm_no);
    DD  = CalculateDDHc(Z, nviews, i_layer);
    SigmaHvd = CalculateSigmaHvd(Hv, nviews, i_layer);

    DXPlus = (abs(DX) + DX)/2;
    DXMinus = (abs(DX) - DX)/2;

    DDPlus = (abs(DD) + DD)/2;
    DDMinus = (abs(DD) - DD)/2;

    if alpha ~= 0 && delta ~=0
        HcUp = DXPlus + DDMinus*(Hc{1,i_layer} + SigmaHvd) + alpha*Hc{1,i_layer}*LcMinus + delta*Hc{1,i_layer};
        HcUn = DXMinus + DDPlus*(Hc{1,i_layer} + SigmaHvd) + alpha*Hc{1,i_layer}*LcPlus;
    elseif alpha ~= 0
        HcUp = DXPlus + DDMinus*(Hc{1,i_layer} + SigmaHvd) + alpha*Hc{1,i_layer}*LcMinus;
        HcUn = DXMinus + DDPlus*(Hc{1,i_layer} + SigmaHvd) + alpha*Hc{1,i_layer}*LcPlus;
    elseif delta ~= 0
        HcUp = DXPlus + DDMinus*(Hc{1,i_layer} + SigmaHvd) + delta*Hc{1,i_layer};
        HcUn = DXMinus + DDPlus*(Hc{1,i_layer} + SigmaHvd);
    else
        HcUp = DXPlus + DDMinus*(Hc{1,i_layer} + SigmaHvd);
        HcUn = DXMinus + DDPlus*(Hc{1,i_layer} + SigmaHvd);
    end

    % Equation 9
    Hc{1,i_layer} = Hc{1,i_layer}.*power((HcUp./HcUn),(0.5));
    Hc{1,i_layer}(isnan(Hc{1,i_layer})|isinf(Hc{1,i_layer}))=0;

end

function Hv = updateHv(R, Z, Hc, Hv, beta, v_ind, i_layer, LvdPlus, LvdMinus)
    
    DX = CalculateThetaTranspose(Z, v_ind, i_layer) * R;
    DD  = CalculateThetaTranspose(Z, v_ind, i_layer) * CalculateTheta(Z, v_ind, i_layer);

    DXPlus = (abs(DX) + DX)/2;
    DXMinus = (abs(DX) - DX)/2;

    DDPlus = (abs(DD) + DD)/2;
    DDMinus = (abs(DD) - DD)/2;
    
    if beta ~= 0
        HvUp = DXPlus + DDMinus*(Hc{1,i_layer} + Hv{v_ind, i_layer}) + beta*Hv{v_ind, i_layer}*LvdMinus{v_ind,1};
        HvUn = DXMinus + DDPlus*(Hc{1,i_layer} + Hv{v_ind, i_layer}) + beta*Hv{v_ind, i_layer}*LvdPlus{v_ind,1};
    else
        HvUp = DXPlus + DDMinus*(Hc{1,i_layer} + Hv{v_ind, i_layer});
        HvUn = DXMinus + DDPlus*(Hc{1,i_layer} + Hv{v_ind, i_layer});
    end

    % Equation 12
    Hv{v_ind, i_layer} = Hv{v_ind, i_layer}.*power((HvUp./HvUn),(0.5));
    Hv{v_ind, i_layer}(isnan(Hv{v_ind, i_layer})|isinf(Hv{v_ind, i_layer}))=0;   
end

function error = cost_function(X, Z, Hc, Hv, Lc, Lvd, nviews, nOfLayers, alpha, beta, delta)
    
    part_one = 0;
    for v_ind = 1:nviews
        R = transpose(X{1,v_ind});
        %R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
        part_one = part_one + norm(R - (CalculateTheta(Z, v_ind, nOfLayers)*(Hc{1,nOfLayers} + Hv{v_ind, nOfLayers})), 'fro');
        
    end
    
    part_two = 0;
    if alpha ~=0 
        for i_layer = 1:nOfLayers

            part_two = part_two + trace(Hc{1,i_layer}*Lc*transpose(Hc{1,i_layer}));
        end
        part_two = alpha*part_two;
    end 
    
    part_three = 0;
    if beta ~=0
        for v_ind = 1:nviews
            for i_layer = 1:nOfLayers

                part_three = part_three + trace(Hv{v_ind, i_layer}*Lvd{v_ind,1}*transpose(Hv{v_ind, i_layer}));
            end
        end
        part_three = beta*part_three;
    end
    
    part_four = 0;
    if delta ~=0 
        for i_layer = 1:nOfLayers

            part_four = part_four + trace(Hc{1,i_layer}*transpose(Hc{1,i_layer}));
        end
        part_four = delta*part_four;
    end 
    
    % Equation 6
    error = part_one + part_two + part_three - part_four;
    
end

function DXZ = CalculateDXZ(Z, R, v_ind, i_layer)
    ThetaTranspose = CalculateThetaTransposeWithMinus(Z, v_ind, i_layer);
    
    DXZ = ThetaTranspose * R;
end

function DDZ = CalculateDDZ(Z, v_ind, i_layer)

    ThetaTransposeWithMinus = CalculateThetaTransposeWithMinus(Z, v_ind, i_layer);
    Theta = CalculateTheta(Z, v_ind, i_layer);
    
    DDZ = ThetaTransposeWithMinus * Theta;
end

function ThetaTransposeWithMinus = CalculateThetaTransposeWithMinus(Z, v_ind, i_layer)
    ThetaTransposeWithMinus = 0;
    
    for i = (i_layer-1):-1:1
        if i == (i_layer-1)
            ThetaTransposeWithMinus = transpose(Z{v_ind, i});
        else
            ThetaTransposeWithMinus = ThetaTransposeWithMinus * transpose(Z{v_ind, i});
        end
    end
end

function Theta = CalculateTheta(Z, v_ind, i_layer)
    Theta = 0;
    
    for i = 1:i_layer
        if i == 1
            Theta = Z{v_ind, i};
        else
            Theta = Theta * Z{v_ind, i};
        end
    end
end

function ThetaTranspose = CalculateThetaTranspose(Z, v_ind, i_layer)
    ThetaTranspose = 0;
    
    for i = i_layer:-1:1
        if i == i_layer
            ThetaTranspose = transpose(Z{v_ind, i});
        else
            ThetaTranspose = ThetaTranspose * transpose(Z{v_ind, i});
        end
    end
end

function DXHc = CalculateDXHc(Z, X, nviews, i_layer, algorithm_no)
    DXHc = 0;
    
    for v_ind = 1:nviews
        if algorithm_no == 1
            R = abs(transpose(X{1,v_ind}));
        elseif algorithm_no == 2
            R = transpose(X{1,v_ind});
        end
        %R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
        
        DXHc = DXHc + (CalculateThetaTranspose(Z, v_ind, i_layer) * R);
    end
end

function DDHc = CalculateDDHc(Z, nviews, i_layer)
    DDHc = 0;
    
    for v_ind = 1:nviews        
        DDHc = DDHc + (CalculateThetaTranspose(Z, v_ind, i_layer) * CalculateTheta(Z, v_ind, i_layer));
    end
end

function SigmaHvd = CalculateSigmaHvd(Hv, nviews, i_layer)
    SigmaHvd = 0;
    
    for v_ind = 1:nviews
        SigmaHvd = SigmaHvd + Hv{v_ind, i_layer};
    end
end