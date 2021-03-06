function newmap = Small_Block_Filter(map, N, small_size)
%--------------------------------------------------------------------------    
% Small_Block_Filter is a filter to removing the small isolated regions.
%
% Implemented by ZhangYu, Image Processing Center, BUAA
%--------------------------------------------------------------------------
% Input:
%         map: input decision map;
%           N: number of the source images here;
%  small_size: the minimum permitted region size existed in the map;
%
% Output:
%      result: return the filtering result
%--------------------------------------------------------------------------

    % define the size of the small region
    P = small_size;
    for ii = 1 : N
        ptmap = (map == ii);

        % Process the Positive image, delete the small patches
        pL = bwlabel(ptmap, 8);

        % Compute the area of each component.
        pS = regionprops(pL, 'Area');
        % Remove small objects.
        ppL = ismember(pL, find([pS.Area] >= P));

        % Process the Negative image, delete the small patches    
        ntmap = ~(ppL > 0);
        nL = bwlabel(ntmap, 8);

        % Compute the area of each component.
        nS = regionprops(nL, 'Area');
        % Remove small objects.    
        npL = ismember(nL, find([nS.Area] >= P));
        tmap = ~(npL > 0);
        result(:,:,ii) = tmap;
    end

    % Find if there are confused pixels
    sCount = sum(result, 3);
    Tag = (sCount == 1);
    newmap = zeros(size(map));
    for ii = 1 : N
       newmap = newmap + ii * result(:,:,ii) .* Tag;
    end
end