function result = Local_Saliency(g, sz)
%--------------------------------------------------------------------------    
% Local_Saliency is a method to compute the weighted gradients 
% in a local window.
% Implemented by ZhangYu, Image Processing Center, BUAA
%--------------------------------------------------------------------------
% Input:
%           g: input image
%          sz: size of the local window
% Output:
%      result: return the computing result
% Example
%       g = magic(5);
%       sz = 3;
%       ls = Local_Saliency(g, sz);
%--------------------------------------------------------------------------

    % Center
    c = floor(sz / 2);

    % Filer operator
    H = ones(sz, sz);
    center = c + 1;
    for ii = 1 : sz
        for jj = 1 : sz
            dis = sqrt( (ii -center) ^ 2 + (jj - center) ^ 2 );
            H(ii, jj) = 1 / (dis + 1);
        end
    end


    % Compute local saliency
    result = imfilter(g, H, 'symmetric');
end