function branchCutMap = BranchCutsJME(residueMap, maxBranchCutLength, maxSearchWinSize)
%
% Creator John Merryman, DTU Space
% Date: 05 Mar. 2018
%
% Usage:
%
% branchCutMap = BranchCutsJME(residueMap, [maxBranchCutLength], [maxSearchWinSize])
%
% Place branch cuts for phase unwrapping according to:
% Goldstein, Zebker, Werner. "Satellite radar interferometry: two
% dimensional phase unwrapping", Radio Science, 1988.
%
% residueMap         : (input)  Phase residue matrix. Values +1 and -1 are expected for 
%                               positive and negative phase residues. 0 elsewhere.
% maxBranchCutLength : (input)  Maximum branch cut length in pixels.
% branchCutMap       : (output) 1 if branchcut, 0 elsewhere.
%
% Example:
%   residueMap = PhaseResiduesJME(phase);
%   branchCutMap = BranchCutsJME(residueMap);
%

    if nargin < 2
        maxBranchCutLength = 200;
    end

    if nargin < 3
        %maxSearchWinSize = max(size(residueMap,1), size(residueMap,2));
        maxSearchWinSize = 200;
    end       
    
    % Inits
    branchCutMap = zeros(size(residueMap));
    
    % Place branch cuts where there are nans in residue map.
    branchCutMap(isnan(residueMap)) = 1;

    % Find residues    
    [residueRows, residueCols, residueCharge] = find(residueMap);
    if isempty(residueRows) || isempty(residueCols)
        fprintf('WARNING: No residues found\n');
        return
    end
    residueRows   = residueRows(~isnan(residueCharge));
    residueCols   = residueCols(~isnan(residueCharge));
    residueCharge = residueCharge(~isnan(residueCharge));    
    
    % For every residue, try to form discharged trees in increasingly large
    % neighbourhoods (search window sizes)
    numOfResidues = size(residueRows,1);
    wasDischarged = zeros(numOfResidues, 1);
    idx = (1 : numOfResidues)';
    listOfResidues = [residueRows residueCols residueCharge wasDischarged idx];

    fprintf('numOfResidues     = %d\n', numOfResidues)

    treeCount = 0;    
    for searchWinSize = 3 : 2 : maxSearchWinSize
        numOfChargedTrees = 0;
        for i = 1 : numOfResidues
            % If not already discharged ...
            if (listOfResidues(i, 4) == 0)
                treeCount = treeCount + 1;
                %fprintf('Growing tree %d: \n', treeCount);
                
                % Try to grow an uncharged tree using this window size.
                % If tree has neutral charge, update listOfResidues marking all the residues in it as discharged.
                [residueTree, listOfResidues] = ...
                    growTree(i, listOfResidues, size(residueMap), searchWinSize, maxBranchCutLength);
                %disp(residueTree);
                
                treeCharge = calcResidueTreeCharge(residueTree);
                if treeCharge ~= 0
                   numOfChargedTrees = numOfChargedTrees + 1;
                   %fprintf('WARNING: Tree %d has non-zero charge: %d\n', treeCount, treeCharge);
                   %disp(residueTree);
                end

                %fprintf('numOfNodes = %d\n', size(residueTree, 1));
                %if size(residueTree, 1) > 2
                %disp(residueTree);    
                %end
                %listOfResidues(residueTreeIndices)

                % Connect nodes in the tree and update branch cut map
                branchCutMap = connectPoints(branchCutMap, residueTree(:,1:2));

            end
        end % end loop over residues
                
        %if (numOfChargedTrees ~= 0)
        %    fprintf('Search win. size: %d, numOfChargedTrees = %d\n', searchWinSize, numOfChargedTrees);    
        %else
        %    fprintf('Search win. size: %d. All trees are discharged.\n', searchWinSize);    
        %    break
        %end
        
    end % end loop over search window sizes
    
    fprintf('numOfResidueTrees = %d\n', treeCount)
    
end

function [residueTree, listOfResidues] = ...
    growTree(currentResidueIdx, listOfResidues, residueMapSize, searchWinSize, maxTreeLength)
%
% Grow an uncharged residue tree according to:
% Goldstein, Zebker, Werner. "Satellite radar interferometry: two
% dimensional phase unwrapping", Radio Science, 1988.
%
% currentResidueIdx   : (input)  Integer index of starting residue
% listOfResidues      : (input)  N x 3 vector with residue row, column and charges
% residueMapSize      : (input)  1 x 2 vector with residue map number of rows and columns
% searchWinSize       : (input)  Search window size in pixels.
% maxTreeLength       : (input)  Maximum branch cut length in pixels.
% residueTreeIndices  : (output) N x 2 list of residues forming the tree (rows and columns).

    if nargin < 4
        searchWinSize = 3;
    end    

    if nargin < 5
        maxTreeLength = 50;
    end
    
    % Init tree with current residue
    residueTree = listOfResidues(currentResidueIdx, :);  
    treeCharge = calcResidueTreeCharge(residueTree);
    
    %For now mark current residue as visited, to exclude it from neighbour
    %searches
    listOfResidues(currentResidueIdx, 4) = 1;
    
    treeLength = 0;
    % While tree is charged ...
    while (treeCharge ~= 0)
        
        % Search for a new tree node.
        % The nodes can be residues or border pixels.
        [newTreeNode, visitedResidueIdx, treeLengthIncrease] = ...
            findNewTreeNode(residueTree, listOfResidues, residueMapSize, searchWinSize);
    
        % If new node was found ...
        if ~isempty(newTreeNode)
            % If within max. tree length ...
            %fprintf('Adding new node ...\n')
            %fprintf('Current tree length: %d\n', treeLength)
            treeLength = treeLength + treeLengthIncrease;
            if (treeLength <= maxTreeLength)
                % Add node to tree 
                residueTree = [residueTree; newTreeNode];
                % Update tree charge
                treeCharge = calcResidueTreeCharge(residueTree);
                % If the node was a residue, mark it as visited for now 
                % to exclude it from further neighbour searches.
                if (visitedResidueIdx ~= 0)
                   listOfResidues(visitedResidueIdx, 4) = 1;
                end
            % If tree would become too long, rather leave it charged and try
            % later with a larger search window.
            else
                break
            end
        % If no new node found, leave tree charged                
        else
            break
        % If not found within max search window, choose closest image 
        % border pixel and consider tree discharged
            %lastResidueCoords = listOfResidues(visitedResidueIndices(end), 1:2);
            %[closestBorderCoords, d] = findClosestImageBorderPixel(lastResidueCoords, residueMapSize);
            %closestBorderCoords = [closestBorderCoords 0 1];
            %residueTree = listOfResidues(visitedResidueIndices, :);
            %residueTree = [residueTree; closestBorderCoords];
        end

    end
    
    % If tree is charged, 
    if (treeCharge ~= 0)
        for i = 1 : size(residueTree, 1)
            % If tree node is a residue ...
            if (residueTree(i,5) ~= 0)
                % Mark all residues as charged, so they can be
                % searched again with larger window sizes.
                listOfResidues(residueTree(i,5), 4) = 0;
            end
        end
    end
    
end

function treeCharge = calcResidueTreeCharge(residueTree)
    treeCharge = sum(residueTree(:,3));
end

function [newTreeNode, visitedResidueIdx, treeLengthIncrease] = ...
    findNewTreeNode(residueTree, listOfResidues, residueMapSize, searchWinSize)

    newTreeNode = [];
    visitedResidueIdx = 0;
    treeLengthIncrease = 0;
    residueSearchResults = [];
        
    % Define search intervals
    lastResidueRow = residueTree(end, 1);
    lastResidueCol = residueTree(end, 2);
    maxRowIndex = residueMapSize(1);
    maxColIndex = residueMapSize(2);
    searchWin2 = floor(searchWinSize/2);

    rowSearchIntervalStart = max(1, lastResidueRow - searchWin2);
    rowSearchIntervalStop  = min(maxRowIndex, lastResidueRow + searchWin2);
    colSearchIntervalStart = max(1, lastResidueCol - searchWin2);
    colSearchIntervalStop  = min(maxColIndex, lastResidueCol + searchWin2);

    treeCharge = calcResidueTreeCharge(residueTree);       
        
    % Find new residue candidates not included in this or previous trees.
    unvisitedIndices = find(listOfResidues(:,4) == 0);
    if ~isempty(unvisitedIndices)
        listOfUnvisitedResidues = listOfResidues(unvisitedIndices, :);
        residueSearchResults = find(listOfUnvisitedResidues(:,1) >= rowSearchIntervalStart & ...
                                    listOfUnvisitedResidues(:,1) <= rowSearchIntervalStop & ...
                                    listOfUnvisitedResidues(:,2) >= colSearchIntervalStart & ...
                                    listOfUnvisitedResidues(:,2) <= colSearchIntervalStop);
    end             
                     
    % 1st choice: residues which can (help) discharge the tree.
    if ~isempty(residueSearchResults)
        % Get indices in original residue list
        nodeCandidateIndices = unvisitedIndices(residueSearchResults);       
               
        nodeCandidateCharges = listOfResidues(nodeCandidateIndices, 3);        
        withOppositeCharge = find(treeCharge * nodeCandidateCharges < 0);
        if ~isempty(withOppositeCharge)
            %fprintf('DEBUG: choice 1\n');
            visitedResidueIdx = nodeCandidateIndices(withOppositeCharge(1));
            listOfResidues(visitedResidueIdx, 4) = 1; % Mark as visited
            newTreeNode = listOfResidues(visitedResidueIdx, :);
        end
    end

    % 2nd choice: border  pixel which can discharge the tree.
    if isempty(newTreeNode)
        %fprintf('DEBUG: choice 2\n');
        [closestBorderCoords, d] = findClosestImageBorderPixel([lastResidueRow lastResidueCol], residueMapSize);
        if (d <= searchWin2)
            newTreeNode = [closestBorderCoords -treeCharge 1 0];
        end
    end      
    
    % 3rd choice: residues which will not (help) discharge tree.
    if isempty(newTreeNode) && ~isempty(residueSearchResults)
        %fprintf('DEBUG: choice 3\n');
        visitedResidueIdx = nodeCandidateIndices(1);
        listOfResidues(visitedResidueIdx, 4) = 1; % Mark as visited
        newTreeNode = listOfResidues(visitedResidueIdx, :);
    end
        
    % Estimate tree length increase
    if ~isempty(newTreeNode)
        newResidueRow = newTreeNode(1);
        newResidueCol = newTreeNode(2);
        treeLengthIncrease = sqrt((newResidueRow - lastResidueRow)^2 + ...
                                  (newResidueCol - lastResidueCol)^2);
    end                 
end

function [closestBorderCoords, closestDistance] = findClosestImageBorderPixel(residueCoords, residueMapSize)

    residueRow = residueCoords(1);
    residueCol = residueCoords(2);
    
    borderPointRows = [residueRow residueRow 1 residueMapSize(1)];
    borderPointCols = [1 residueMapSize(2) residueCol residueCol];

    d = sqrt((borderPointRows - residueRow).^2 + (borderPointCols - residueCol).^2);
    [M, closestIdx] = min(d);
    closestBorderCoords = [borderPointRows(closestIdx) borderPointCols(closestIdx)];
    closestDistance = d(closestIdx);
end

function inputImage = connectPoints(inputImage, listOfPoints)
    
    if size(listOfPoints,1) < 1
        error('At least one points required.\n')
    end
    
    if size(listOfPoints,2) ~=  2
        error('List of points must include x and y coords.\n')
    end
        
    % Connect each point in the list with the next
    numOfPoints = size(listOfPoints,1);
    
    if numOfPoints == 1
        inputImage(listOfPoints(1,1),listOfPoints(1,2)) = 1;
    else
        for p = 1 : numOfPoints - 1
            row1 = listOfPoints(p,1);
            row2 = listOfPoints(p+1,1);
            col1 = listOfPoints(p,2);
            col2 = listOfPoints(p+1,2);
            d = round(sqrt((row2 - row1)^2 + (col2 - col1)^2));
            for n = 0 : 1/d : 1 
                rown = round(row1 +(row2 - row1) * n); 
                coln = round(col1 +(col2 - col1) * n); 
                inputImage(rown,coln) = 1; 
            end
        end
    end
    
end