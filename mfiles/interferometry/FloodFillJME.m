function unwrappedPhase = FloodFillJME(wrappedPhase, branchCuts, refRowCol)
%
% Creator John Merryman, DTU Space
% Date: 05 Mar. 2018
%
% Usage:
%
% unwrappedPhase = FloodFillJME(wrappedPhase, branchCuts, refRowCol)
%
% Integrate wrapped phase spatially avoiding to cross branch cuts.
% NOTE: if phase matrix size is nrows x ncols -> branchCuts is expected to be 
%       nrows + 1 x ncols + 1
%
% phase      : (input)  Wrapped interferometric phase matrix (real values in rad)where.
% branchCuts : (input)  Map of branch cuts (1 where branch cuts. NaNs are considered branch cuts as well).
% refRowCol  : (input)  [startingRow startingColumn]. Phase are unwrapped
%                       with respect to this pixel.
%
% Example:
%    pha = angle(I);
%    residues = PhaseResiduesJME(pha);     
%    branchCuts = BranchCutsJME(residues);
%    phu = FloodFillJME(pha, branchCuts, [160 200]);  
%
    if nargin < 2
      error('Insufficient number of input parameters, exiting.');
    end

    if nargin < 3
      refRowCol = [floor(size(branchCuts,1)/2) floor(size(branchCuts,2)/2)];
    end

    [nrows, ncols] = size(wrappedPhase);

    unwrappedPhase = nan(nrows, ncols);
    
    if isnan(wrappedPhase(refRowCol(1),refRowCol(2)))
        error('Unwrapping seed has a NaN value');
    else
        wrappedPhase = wrappedPhase - wrappedPhase(refRowCol(1),refRowCol(2));
    end

    newSeedList = refRowCol;
    unwrappedPhase(newSeedList(1,1), newSeedList(1,2)) = ...
        wrappedPhase(newSeedList(1,1), newSeedList(1,2));  
    isUnwrapped = zeros(nrows, ncols);
    isUnwrapped(newSeedList(1,1), newSeedList(1,2)) = 1;   
    
    while ~isempty(newSeedList)
        currentSeed = [newSeedList(1,1) newSeedList(1,2)];
        %fprintf('Current seed: (%d,%d)\n', currentSeed(1), currentSeed(2));     
        seedNeighbourRows = [currentSeed(1) - 1, currentSeed(1), currentSeed(1), currentSeed(1) + 1];
        seedNeighbourCols = [currentSeed(2), currentSeed(2) - 1, currentSeed(2) + 1, currentSeed(2)];
        % For every seed neighbour (4-point connectivity)
        for i = 1 : 4
            %fprintf('Considering neighbour: (%d,%d)\n', seedNeighbourRows(i), seedNeighbourCols(i));  
            neighbourSeed = [seedNeighbourRows(i), seedNeighbourCols(i)];
            % If the neighbour was not already unwrapped and it is not a NaN
            if isValidImageIdx(neighbourSeed, [nrows, ncols]) && ...
               ~isUnwrapped(neighbourSeed(1), neighbourSeed(2)) && ...
               ~isnan(wrappedPhase(neighbourSeed(1), neighbourSeed(2)))

                %fprintf('Neighbour has valid phase value.\n');
          
                wrappedPhaseNeighbour = wrappedPhase(neighbourSeed(1), neighbourSeed(2));
           
                % Unwrap phase diff w.r.t. seed unless there is a branch-cut in between
                % seed and neighbour
                isBranchCut = branchCutInBetween(currentSeed, neighbourSeed, branchCuts);
                
                if (isBranchCut == 0)
                    %fprintf('Unwrapping neighbour ...\n');  
                    phu = unwrap([unwrappedPhase(currentSeed(1), currentSeed(2)) wrappedPhaseNeighbour]');
                    unwrappedPhase(neighbourSeed(1), neighbourSeed(2)) = phu(2);
                    %fprintf('Wrapped phase (seed, neighbour)   = (%f,%f) cycles\n', wrappedPhaseSeed / 2 / pi, wrappedPhaseNeighbour / 2 / pi);
                    %fprintf('Unwrapped phase (seed, neighbour) = (%f,%f) cycles\n', unwrappedPhase(currentSeed(1), currentSeed(2)) / 2 / pi, unwrappedPhase(neighbourSeed(1), neighbourSeed(2)) / 2 / pi);
                    % Mark as unwrapped
                    isUnwrapped(neighbourSeed(1), neighbourSeed(2)) = 1;
                    % add it to candidate seed stack
                    newSeedList(end + 1,:) = neighbourSeed;
                end
                    
            end
        end
        % Remove current seed from stack
        newSeedList(1,:) = [];
    end

end

function isValid = isValidImageIdx(rowCol, imageSize)
    isValid = 0;
    if rowCol(1) >= 1 && rowCol(1) <= imageSize(1) && ...
       rowCol(2) >= 1 && rowCol(2) <= imageSize(2)
        isValid = 1;
    end
end

function isBranchCut = branchCutInBetween(rowCol1, rowCol2, branchCuts)
    isBranchCut = 0;
    % If the two pixels are not the same
    if ~isempty(find(diff([rowCol1; rowCol2]), 1))
        % If pixels are on the same row
        if rowCol1(1) == rowCol2(1)
            col1 = rowCol1(2);
            col2 = rowCol2(2);
            row  = rowCol1(1);
            % Check for vertical branch cuts
%             if (branchCuts(row, max(col1,col2))         == 1 || ...
%                 branchCuts(row + 1, max(col1,col2))     == 1) ...
%                 isBranchCut = 1;
%             end
            
            if (branchCuts(row, max(col1,col2))         == 1 || ...
                branchCuts(row + 1, max(col1,col2))     == 1 || ...
                branchCuts(row - 1, max(col1,col2))     == 1 || ...
                branchCuts(row + 2, max(col1,col2))     == 1 || ...
                branchCuts(row, max(col1,col2) - 1)     == 1 || ...
                branchCuts(row + 1, max(col1,col2) - 1) == 1 || ...
                branchCuts(row - 1, max(col1,col2) - 1) == 1 || ...
                branchCuts(row + 2, max(col1,col2) - 1) == 1 || ...
                branchCuts(row, max(col1,col2) + 1)     == 1 || ...
                branchCuts(row + 1, max(col1,col2) + 1) == 1 || ...
                branchCuts(row - 1, max(col1,col2) + 1) == 1 || ...
                branchCuts(row + 2, max(col1,col2) + 1) == 1)
                isBranchCut = 1; 
            end
            % Check diagonal branch cuts
            %if (branchCuts(row, max(col1,col2) + 1) == 1 && branchCuts(row + 1, max(col1,col2))     == 1) || ...
            %   (branchCuts(row, max(col1,col2))     == 1 && branchCuts(row + 1, max(col1,col2) + 1) == 1)
            %    isBranchCut = 1; 
            %end
        % If pixels are on the same column
        elseif rowCol1(2) == rowCol2(2)
            % Check horizontal branch cuts
            row1 = rowCol1(1);
            row2 = rowCol2(1);
            col  = rowCol1(2);
%             if (branchCuts(max(row1,row2), col)     == 1 || ...
%                 branchCuts(max(row1,row2), col + 1) == 1)
%                 isBranchCut = 1; 
%             end
            if (branchCuts(max(row1,row2), col)     == 1 || ...
                branchCuts(max(row1,row2), col + 1) == 1 || ...
                branchCuts(max(row1,row2), col - 1) == 1 || ...
                branchCuts(max(row1,row2), col + 2) == 1 || ...
                branchCuts(max(row1,row2) - 1, col)     == 1 || ...
                branchCuts(max(row1,row2) - 1, col + 1) == 1 || ...
                branchCuts(max(row1,row2) - 1, col - 1) == 1 || ...
                branchCuts(max(row1,row2) - 1, col + 2) == 1 || ...
                branchCuts(max(row1,row2) + 1, col)     == 1 || ...
                branchCuts(max(row1,row2) + 1, col + 1) == 1 || ...
                branchCuts(max(row1,row2) + 1, col - 1) == 1 || ...
                branchCuts(max(row1,row2) + 1, col + 2) == 1)
                isBranchCut = 1; 
            end
            % Check diagonal branch cuts
            %if (branchCuts(max(row1,row2) + 1, col) == 1 && branchCuts(max(row1,row2), col + 1)     == 1) || ...
            %   (branchCuts(max(row1,row2), col)     == 1 && branchCuts(max(row1,row2) + 1, col + 1) == 1)
            %    isBranchCut = 1; 
            %end
        end
    end
end