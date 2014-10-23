function result = block_values(blocks, S, maxDim, value)
dim = maxDim;
while dim > 0    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)
    [values, Sind] = qtgetblk(blocks, S, dim);
    values(end,1:end,:) = value;
    values(1:end,end,:) = value;
    values(1,1:end,:) = value;
    values(1:end,1,:) = value;
    blocks = qtsetblk(blocks,S,dim,values);
  end
  dim = dim / 2;
end

blocks(end,1:end) = value;
blocks(1:end,end) = value;
blocks(1,1:end) = value;
blocks(1:end,1) = value;

result = blocks;