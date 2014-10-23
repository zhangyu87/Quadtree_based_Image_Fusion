function [S, Fusion_tag, maxFM] = qt_sm_fusion(Imgs, FMs, Num, Def_level)
%------------------------------------------------------------------------
% qt_sm_fusion is the function to decompose the input images 
% in a quadtree structure
% 


[M,N] = size(Imgs(:,:,1));
% maxDim = M / 2;
maxDim = M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing the minFM, maxFM and sumFM at each pixel for the whole map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minFM = FMs(:,:,1);
maxFM = minFM;
for ii = 2 : Num
minFM = min(minFM, FMs(:,:,ii));
maxFM = max(maxFM, FMs(:,:,ii));
end

%% Initialize blocks
S = zeros(M,N);
S(1:maxDim:M, 1:maxDim:N) = maxDim;

%% Quad-tree decompose process
level = 1;
dim = maxDim;
maxlevel = floor(log2(M));
maxlevel = min(maxlevel, Def_level);

%% Decomposition of the source images
while (level < maxlevel)	
    [blockVal0, Sind] = qtgetblk(FMs(:,:,1), S, dim);
	%% End condition
    if (isempty(Sind))
        % Done!
        break;
    end
    %% Split process
    doSplit = QTDECOMP_Split(S, FMs, dim, Num, minFM, maxFM);
    %% Next level
    dim = dim/2;
	%% Split blocks for doSplit == 1
    Sind = Sind(doSplit);
    %% Split the blocks
    Sind = [Sind ; Sind+dim ; (Sind+M*dim) ; (Sind+(M+1)*dim)];
    S(Sind) = dim;
    %% Decomposit in the next level
    level = level + 1;
end

%% Initialise data
S = sparse(S);
dim = maxDim;
level = 1;
%% Define a tag for most FM images
FM_tag = zeros(M,N);

%% Begin loop
while (level <= maxlevel)

    len = length(find(S == dim));
    
    if len ~= 0
        
        values = zeros(Num, len);
        for i = 1 : Num
            [blockFM, Sind] = qtgetblk(FMs(:,:,i), S, dim);
        %   focus-measure of the blocks
            values(i,:) = sum(sum(blockFM));
        end
        % calculate the maximum focus-measure of the corresponding blocks
        maxVals = max(values);
        % Find the index of the blocks with maximum focus-measure
        is_equal = zeros(1, len);
        eq_count = is_equal;
        index = is_equal;

        for i = 1 : Num
           is_equal = (maxVals == values(i, : ));
           eq_count = eq_count + is_equal;
           index = index + is_equal * i;
        end
        %% Find the image index
        index = index .* (eq_count == 1) + (-1) * (eq_count > 1);

        %% Judge which image has the greatest FM
        blocktag = zeros(size(blockFM));
        for i = 1 : len
            blocktag( : , : , i) = index(1, i);
        end
        %% Set back the fusion tag
        FM_tag = qtsetblk(FM_tag,S,dim,blocktag);
    end
	%% Next Scale
	dim = dim/2;
    level = level + 1;
end	

Fusion_tag = FM_tag;