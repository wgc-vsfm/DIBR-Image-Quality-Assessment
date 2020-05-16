load DIBR\DIBR_dmos
% tic
% % DIBR_apt = our_DIBR(imgB);
% % toc
cnt = 0;
num = 84;
DIBR_apt = zeros(num,1);
for cnt=1:num
cnt
img1_name = [DIBR_name{cnt,1}(19:end)];
img2_name = [DIBR_name{cnt,2}(19:end)];
imgA = imread(img1_name);
imgB = imread(img2_name); 
DIBR_apt(cnt) = our_DIBR(imgB);
end
% toc
% load A3_MOS.mat
% load A3_OUR.mat
% DIBR_apt=A1_OUR;
% DIBR_dmos=A1_MOS;
% % % load SSIM_DIBR
% % DIBR_apt=dibr_biqme;
cc = abs(corr(DIBR_apt,DIBR_dmos,'type','spearman'));
[delta,beta,x,y,diff] = findrmse2(DIBR_apt,DIBR_dmos);
[corr(x,y),corr(x,y,'type','spearman'),corr(x,y,'type','kendall'),mean(abs(diff)),(mean(diff.^2))^0.5]
