function result = SML(g, step, T)
%--------------------------------------------------------------------------    
% SML is the sum of modified Laplacian, which is well-known focus measure 
% for multifocus image fusion.   
% Implemented by ZhangYu, Image Processing Center, BUAA
%--------------------------------------------------------------------------
% Input:
%       g: input matrix
%    step: sample distance 
%       T: thresholding value
% Output:
%  result: return the 
% Example
%       g = magic(10);
%    step = 1;
%       T = 5;
%  result = = SML(g, step, T);
%--------------------------------------------------------------------------

    % padding the source matrix for convenient computation
    domain = ones(2 * step);
    gpad = hPadImage(g, domain, 'symmetric');

    [M, N] = size(g);
    result = zeros(M, N);

    % A basic way to compute SML
% %     for i = 1 : M 
% %         for j = 1 : N
% %             x = i + step;
% %             y = j + step;
% %             part1 = abs(2 * gpad(x, y) - gpad(x - step, y) - gpad(x + step, y));
% %             part2 = abs(2 * gpad(x, y) - gpad(x, y - step) - gpad(x, y + step));
% %             V = part1 + part2;
% %             
% %             result(i,j) = V;
% %         end
% %     end
    
    % A fast way to compute the SML by using the matrix operations
    mat1 = gpad(step + 1 : M + step, step + 1 : N + step);
    
    mat2 = gpad(1 : M, step + 1 : N + step);
    mat3 = gpad(2 * step + 1 : M + 2 * step, step + 1 : N + step);
    
    mat4 = gpad(step + 1 : M + step, 1 : N);
    mat5 = gpad(step + 1 : M + step, 2 * step + 1 : N + 2 * step);
    
    % Two components
    part1 = abs(2 * mat1 - mat2 - mat3);
    part2 = abs(2 * mat1 - mat4 - mat5);
    
    result = part1 + part2;  
    
    % Threshol
    result(result<T) = 0;
end