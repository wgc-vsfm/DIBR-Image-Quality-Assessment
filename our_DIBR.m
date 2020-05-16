function score=our_DIBR(org_img)
if (size(org_img,3)==3),
    org_img = double(rgb2gray(org_img));
else org_img = double(org_img);
end
tt=1;
%小波分解
[cA,cH,cV,cD]=dwt2(org_img,'db20');
%提取blur
ss_cA=ssq(cA);
ss_cH=ssq(cH);
ss_cV=ssq(cV);
ss_cD=ssq(cD);
%ss=ss_cD +(ss_cH+ss_cV)/2 + 0.9*ss_cA;
%ss=ss_cD;
ss=0.5*ss_cD +0.3*(ss_cH+ss_cV)/2+0.2*ss_cA;
%dst = max(sum(ss));
%提取几何失真
thresh = graythresh(cA);     %自动确定二值化阈值
I2_cA = im2bw(cA,thresh); 
BW1_cA = edge(I2_cA,'canny'); 
BW1_cH = edge(cH,'canny');
BW1_cV = edge(cV,'canny');
BW1_cD = edge(cD,'canny');
G1_cH = 2*BW1_cA.*BW1_cH;
G2_cH = BW1_cA.^2 + BW1_cH.^2;
yy_cH = (G1_cH + tt) ./(G2_cH + tt);
yy_cH=yy_cH./max(yy_cH(:));
G1_cV = 2*BW1_cA.*BW1_cV;
G2_cV = BW1_cA.^2 + BW1_cV.^2;
yy_cV = (G1_cV + tt) ./(G2_cV + tt);
yy_cV=yy_cV./max(yy_cV(:));
G1_cD = 2*BW1_cA.*BW1_cD;
G2_cD = BW1_cA.^2 + BW1_cD.^2;
yy_cD = (G1_cD + tt) ./(G2_cD + tt);
yy_cD=yy_cD./max(yy_cD(:));
jihe = sum(sum(yy_cH))+sum(sum(yy_cV))+sum(sum(yy_cD));
% %图像复杂度
org_img = imresize(org_img,1/4);
fe = complexity_index(org_img);
%质量分数
%score=ss;
%score= (0.3*(ss/2.4185) + jihe/616577)/fe;
score1=(ss/2.4185)/fe;
score2=(jihe/616577)/fe;
a=0.15;
score=(1/(1+a))*score2+(a/(1+a))*score1;
function ss = ssq(bands)
lh_img = bands.^2;
ss = log10(1+mean(lh_img(:)));
function fe = complexity_index(im)
ar = armodel(im,1);
bi = bfilter(im/255,1,[3,0.1],1)*255;
fe = fe_index(0.1*ar+0.9*bi,im,1);
function ff = fe_index(img1,img2,tt)
ee = img2-img1;
xx = -255:255;
yy = round(ee(1:tt:end,1:tt:end));
[nn,~] = hist(yy(:),xx);
pp = (1+2*nn)/sum(1+2*nn);
ff = -sum(pp.*log2(pp));
%=======================================================
function B = bfilter(A,w,sigma,tt)
sigma_d = sigma(1);
sigma_r = sigma(2);
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));
dim = size(A);
B = zeros(dim);
for i = 1:tt:dim(1)
for j = 1:tt:dim(2)
iMin = max(i-w,1);
iMax = min(i+w,dim(1));
jMin = max(j-w,1);
jMax = min(j+w,dim(2));
I = A(iMin:iMax,jMin:jMax);
H = exp(-(I-A(i,j)).^2/(2*sigma_r^2));
F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
B(i,j) = sum(F(:).*I(:))/sum(F(:));
end
end
%=======================================================
function imgrec = armodel(imgin,tt)
sr=3;%search range
mr=1;%model range
imgt=padarray(imgin,[sr+mr sr+mr],'symmetric');
imgrec=zeros(size(imgin));
[m n]=size(imgt);
N=(2*sr+1)^2-1;
K=(2*mr+1)^2-1;
A=zeros(N,K+1);
for ii=mr+sr+1:tt:m-sr-mr
for jj=mr+sr+1:tt:n-sr-mr
con=1;
patch0=imgt(ii-mr:ii+mr,jj-mr:jj+mr);
for iii=-sr:+sr
for jjj=-sr:+sr
if iii==0&&jjj==0
continue;
end
patch=imgt(ii+iii-mr:ii+iii+mr,jj+jjj-mr:jj+jjj+mr);
vec=patch(:);
A(con,:)=vec';
con=con+1;
end
end
b=A(:,mr*(2*mr+2)+1);
A2=A;
A2(:,mr*(2*mr+2)+1)=[];
if rcond(A2'*A2)<1e-7
a = ones(K,1)/K;
else
a = A2\b;
end
vec0=patch0(:);
vec0(mr*(2*mr+2)+1)=[];
rec=vec0'*a;
imgrec(ii-sr-mr,jj-sr-mr)=rec;
end
end