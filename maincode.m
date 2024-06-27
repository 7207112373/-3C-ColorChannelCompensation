% disp('Color Channel Compensation (3C): A Fundamental
% Pre-Processing Step for Image Enhancement')
clc;
close all;
clear all;
warning off
%%%%%under water image enhncement

%%%%%%%%%%%%%%%% Image Aquisition %%%%%%%%%%%%%%%%%%%%%%%%%%%
[filename, pathname] = uigetfile({'*.jpg';'*.png'}, 'pick an image');
if isequal(filename, 0) || isequal(pathname, 0)   
    helpdlg('Image input canceled.');  

else
    X=imread(fullfile(pathname, filename));
end

figure,imshow(X);title('original image');
tic
%%%%%%%%%%%%devide color space images CCC process
  N = 256;
    A = im2uint8(X);
    fig = figure;
    %[peaksnr]= measerr(A,fig);
    %fprintf('The Peak-SNR value is %.4f\n', peaksnr);

    subplot(2,2,1);
    imshow(A);
    title('original');
    ColorList = { 'Red' 'Green' 'Blue' };
    gr = 0:1/(N-1):1;        
           
    for k = 1:3
        
        % color map:
        cMap = zeros(N,3);
        cMap(:,k) = gr;   
        
        % Display monochromatic image:
        subplot(2,2,k+1);
        imshow(ind2rgb(A(:,:,k),cMap));
        title(ColorList{k});
        
    end

%%% proposed work Color Channel Compensation
tic
normal_thr_limit=0.5;
low_limt=0.002;
up_limit=0.999;
%----------------------------------------------------------------------
under_image=X;
[CONTRAST CCC Correction]=size(X);
%----------------------------------------------------------------------
if Correction==3
    inc_pixel_limit=0.04;dec_pixel_limt=-0.04;
    max_chromatic=rgb2ntsc(under_image); %  grouping pixels data
    mean_adjustment=inc_pixel_limit-mean(mean(max_chromatic(:,:,2)));
    max_chromatic(:,:,2)=max_chromatic(:,:,2)+mean_adjustment*(0.989-max_chromatic(:,:,2));
    mean_adjustment=dec_pixel_limt-mean(mean(max_chromatic(:,:,3)));
    max_chromatic(:,:,3)=max_chromatic(:,:,3)+mean_adjustment*(0.898-max_chromatic(:,:,3));
else
    max_chromatic=double(under_image)./255;
end
%----------------------------------------------------------------------

mean_adjustment=normal_thr_limit-mean(mean(max_chromatic(:,:,1)));
max_chromatic(:,:,1)=max_chromatic(:,:,1)+mean_adjustment*(0.958-max_chromatic(:,:,1));
if Correction==3
    max_chromatic=ntsc2rgb(max_chromatic);
end
%----------------------------------------------------------------------
under_image=max_chromatic.*255;
%--------------------caliculate the max to min pixels----------------------
for k=1:Correction
    arr=sort(reshape(under_image(:,:,k),CONTRAST*CCC,1));
    saliency_min(k)=arr(ceil(low_limt*CONTRAST*CCC));
    luminance_max(k)=arr(ceil(up_limit*CONTRAST*CCC));
end
%----------------------------------------------------------------------
if Correction==3
    saliency_min=rgb2ntsc(saliency_min);
    luminance_max=rgb2ntsc(luminance_max);
end

%----------------------------------------------------------------------

reconstucted_image=(under_image-saliency_min(1))/(luminance_max(1)-saliency_min(1));
figure,imshow(mat2gray(under_image));title('underexposed image enhancement');
    
% % %%%%%%%%%%%%%%%temporal correlation pixels separation

            red_color_correlation = adapthisteq(reconstucted_image(:,:,1));
            green_adapthisteq_green = adapthisteq(reconstucted_image(:,:,2));
            blue_adapthisteq_blue = adapthisteq(reconstucted_image(:,:,3));
            fusion_enhanced = cat(3,red_color_correlation,green_adapthisteq_green,blue_adapthisteq_blue);
            figure,imshow(mat2gray(fusion_enhanced));title('proposed enhanced output image');
            toc
     [ PSNR , MSE , MAXERR , L2RAT ] = measerr( reconstucted_image , fusion_enhanced ) ;
disp('PSNR value ')
disp(PSNR)
            
