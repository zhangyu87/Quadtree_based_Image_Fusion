% Matlab demo code for "Quadtree-based multi-focus image fusion using a weighted focus-measure, Information Fusion 22 (2015) 105¨C118". 
% Implemented by ZhangYu, Image Processing Center, BUAA
% Email: uzeful@163.com

% -------------------------------------------------------------------------
% Clear history and memory
% -------------------------------------------------------------------------
clc;clear all;close all;

% -------------------------------------------------------------------------
% Const String: Next line
% -------------------------------------------------------------------------
NextLine = sprintf('\n');

% -------------------------------------------------------------------------
% Algorithm Start Here: Fusion each set of the multifocus source images 
%
% Image Sets including:
% 1.pepsi
% 2.clock
% 3.lab
% 4.toy
% 5.OpenGL
% 6.disk
% 7. toy
% -------------------------------------------------------------------------

for name = 1 : 7

    % Initialise the maximum split level
    level = 0;
    % Initialise the maximum block size permitted
    pmax = 0;
    % default image number and type
	N = 2; type = '.bmp';
	
    % Set the Image name 
    switch name 
        case 1 
            str = 'clock'; 
        case 2 
            str = 'lab';
        case  3 
            str = 'pepsi'; 
        case 4 
            str = 'OpenGL';
        case 5 
            str = 'flower';
        case 6 
            str = 'disk';
        case 7
            str = 'toy'; N = 3;
    end
    fprintf('Processing Image %s :\n',  str);
    
    % ---------------------------------------------------------------------
    % Choose the focus measure for image fusion

    % opt = 1: using SML as the focus-measure;
    % opt = 2: using SWML as the focus-measure;
    % ---------------------------------------------------------------------

    for opt = 2
        % Get dimensions of multi-focus images
        inImg = double(imread (strcat('example images\',str,num2str(1),type)));
        [p1,p2,p3] = size(inImg);

        % Extend the image dimensions to 2 ^ X;  X <= 512
        dim = max(p1, p2);
        maxDim = 1024;
        while(dim < maxDim || dim == maxDim)
            maxDim = maxDim / 2;
        end
        maxDim = maxDim * 2;

        % if the maxDim == 2048, exit the program
        if (maxDim == 2048)
            disp('Image size is too big, it cannot process!');
            return;
        end

        % Store the multi-focus source images
        mImg=zeros(p1,p2,N);
        for q=1:N
            % Read the image data
            inImg = imread (strcat('example images\',str,num2str(q),type));
            [p1,p2,p3] = size(inImg);

            % Convert image data from uint8 to double
            inImg = double(inImg);

            % Store the temp image
            mImg(:,:,q)=inImg;
        end
        
        % -----------------------------------------------------------------
        % Count the running time
        % -----------------------------------------------------------------
        tic
        
        % -----------------------------------------------------------------
        % Compute the gradient map and Focus measures,
        % opt = 1: SML; 
        % opt = 2: SWML;
        % -----------------------------------------------------------------

        % Weighted Block Size
        sz = 17;
        % Initialize the focus measure maps
        FMs = zeros(p1, p2, N);
        if opt == 1
            disp('Fused method 1: Sum of Modified laplacian');
            for kk = 1 : N
                f = SML(mImg( : , : , kk), 1, 5);
                FMs(:,:,kk) = f;
            end
            method = 'SML';
        elseif opt == 2
            disp('Fused method 2:  Sum of Weighted Modified laplacian');
            for kk = 1 : N
                f = SML(mImg( : , : , kk), 1, 5);
                FM = Local_Saliency(f, sz);
                FMs(:,:,kk) = FM;
            end
            method = 'SWML';
        end

        %------------------------------------------------------------------
        % Quad-tree decomposition
        %------------------------------------------------------------------
        
        %% Extend Imgs and Salmaps to maxDim( 2 ^ X )
        Imgs = zeros(maxDim, maxDim, N);

        % If we put the images on the leftup coner of the extended image;
        % Imgs(1 : p1, 1 : p2, : ) = mImg;
        % Salmaps = zeros(maxDim, maxDim, N);
        % Salmaps(1 : p1, 1 : p2, : ) = FMs;

        % Elseif we center the images and the gradient maps on the extended
        % images or maps
        dx = ceil((maxDim - p1) / 2);
        dy = ceil((maxDim - p2) / 2);

        Imgs(dx + 1 : dx + p1, dy + 1 : dy + p2, : ) = mImg;

        Salmaps = zeros(maxDim, maxDim, N);
        Salmaps(dx + 1 : dx + p1, dy + 1 : dy + p2, : ) = FMs;

        %% Decomposition level and  max permitted block-size must be restricted!!!
        pmaxDim = maxDim / 2;

        % Set the Default level, when not set level
        if level == 0
            level = log2(maxDim);
        end

        % Set the default maximum block size permitted
        if pmax == 0
            pmax = maxDim;
        end

        %% Quad tree docomposition process
        pmin = 2;
        Num = N;
        [S, Fusion_tag, maxFM] = qt_sm_fusion(Imgs, Salmaps, Num, level);

        %------------------------------------------------------------------
        % Reconstruction Part
        %------------------------------------------------------------------

        % if we put the image on the left up coner of the extended image;
        % Fusion_decision_map = Fusion_tag(1 : p1, 1 : p2); 

        % elseif we put the source image on the center of the extended image, 
        % extract the fusion decision map from the extended map
        Fusion_decision_map = Fusion_tag(dx + 1 : dx + p1, dy + 1 : dy + p2);
        maxFM = maxFM(dx + 1 : dx + p1, dy + 1 : dy + p2);

        %% 1. First Filter: Open and Close morphilogical filtering
        Iter = 1;
        Fusion_decision_map = morph_filter(Fusion_decision_map, N, Iter);

        %% 2. Second Filter: Fiter small blocks inside
        smallsz = p1 * p2 / 40;
        Fusion_decision_map = Small_Block_Filter(Fusion_decision_map, N, smallsz);


        %------------------------------------------------------------------
        % Fusion Rule
        %------------------------------------------------------------------        
        
        % Initialize the fusion image FImg
        FImg = zeros(p1,p2);

        %% 1. the defined part, copied directly according to the decision map
        for ii = 1 : Num	
            FImg = FImg + mImg(:,:,ii) .* (Fusion_decision_map == ii);
        end


        %% 2. the non-defined part, copied by maximum selection method
        max_tag = zeros(p1, p2, N);
        img_tag = zeros(p1, p2);
        % Find the maximum FM in each FM
        for ii = 1 : N
            tag = (FMs(:,:,ii) == maxFM);
            max_tag(:,:,ii) = tag;
            img_tag = img_tag + tag .* ii;
        end

        % The nonpart images and maximum selection
        non_part = (Fusion_decision_map < 1);
        nonImgs = mImg;
        part1 = zeros(p1,p2);
        for ii = 1 : N
            nonImgs(:,:,ii) = nonImgs(:,:,ii) .* non_part;
            part1 = part1 + nonImgs(:,:,ii) .* max_tag(:,:,ii);
        end

        % The positions where more than one FM(i) have the maxFM
        max_num = sum(max_tag, 3);

        % The sigle and multiple positions
        single_num = (max_num == 1);
        multi_num = 1 - single_num;

        % If there are more than one FMi equal to maxFM
        part2 = sum(nonImgs, 3) ./ N;

        % As to the whole nonpart
        nonPart = part1 .* single_num + part2 .* multi_num;

        % Final Fused image FImg
        FImg = FImg + nonPart;
        FImg = uint8(FImg);

        
        %------------------------------------------------------------------
        % End Counting time
        %------------------------------------------------------------------
        toc
        t = toc;
        
        % Show the fusion image and the corresponding decision map
        figure,imshow(FImg);
        figure,imshow(mat2gray(Fusion_decision_map))

        %--------------------------------------------------------------------------
        % Store the Fusion decision map and the fusion image
        %--------------------------------------------------------------------------
        imwrite((Fusion_decision_map ./ N), strcat('fusion decision maps\', method, '_', str, '.bmp'));
        imwrite(FImg, strcat('fusion results\', method, '_', str, '.bmp'));
    end
    
    % Press 'Enter' key  to process the next set of images
    disp('Press ''Enter'' key  to process the next set of images');
    pause();
    clc, clear all, close all
%% name 1 : 6 loop
end
