%%% QTDECOMP_Split %%%
function dosplit = QTDECOMP_Split(S, FMs, dim, Num, minFM, maxFM)
%% Get the block datas
[blockV, Sind] = qtgetblk(FMs(:,:,1), S, dim);
sumBlock = sum(sum(blockV));
minBlock = sumBlock;
maxBlock = sumBlock;
for i = 2 : Num
[blockV, Sind] = qtgetblk(FMs(:,:,i), S, dim);
tempBlock = sum(sum(blockV));
minBlock = min(minBlock, tempBlock);
maxBlock = max(maxBlock, tempBlock);
end

%% MDFM: maximum difference in the focus measures
MDFM = abs(minBlock - maxBlock);

%% FM for the pixel by pixel
[minblockV, Sind] = qtgetblk(minFM, S, dim);
[maxblockV, Sind] = qtgetblk(maxFM, S, dim);

maxSum = sum(sum(maxblockV));
minSum = sum(sum(minblockV));
%% SMDG: sum of maximum difference of the gradients
SMDG = abs(maxSum - minSum);
%% Compute the threshold for the block to split
T = SMDG * 0.98;
%% Split
dosplit = (MDFM > 0) &(MDFM < T);