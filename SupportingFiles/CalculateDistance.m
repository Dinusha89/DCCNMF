function [distanceMatrix] = CalculateDistance(sourceMatrix, targetMatrix, distanceType)

    % row represents sample and column represents features
    if(distanceType == "Samples")
        noOfRowsInSourceMatrix = size(sourceMatrix,1);      
        noOfRowsInTargetMatrix = size(targetMatrix,1);
  
        distanceMatrix = zeros(noOfRowsInSourceMatrix,noOfRowsInTargetMatrix);

        for sourceRow = 1:noOfRowsInSourceMatrix
            for targetRow = 1:noOfRowsInTargetMatrix
                distanceMatrix(sourceRow, targetRow) = sqrt(sum((sourceMatrix(sourceRow,:) - targetMatrix(targetRow,:)) .^ 2));
            end   
        end
    
    elseif(distanceType == "Features")
        noOfColumnsInSourceMatrix = size(sourceMatrix,2);
        noOfColumnsInTargetMatrix = size(targetMatrix,2);
        
        distanceMatrix = zeros(noOfColumnsInSourceMatrix,noOfColumnsInTargetMatrix);
        
        for sourceColumn = 1:noOfColumnsInSourceMatrix
            for targetColumn = 1:noOfColumnsInTargetMatrix
                distanceMatrix(sourceColumn, targetColumn) = sqrt(sum((sourceMatrix(:,sourceColumn) - targetMatrix(:,targetColumn)) .^ 2));
            end   
        end
    end

end