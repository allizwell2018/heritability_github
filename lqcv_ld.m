%%
function result = lqcv_ld(U,N,n,sigmaE,p,x)

%lambda = [0.1 1 5 10 15 20 25 30 35 40 45 50 60 70 100];
%lambda=10:5:50;
%lambda = 1e-5:1e-5:1e-4;
MaxIterNum = 1e5;
u_Init = randn(N,1);

% perm_index1=randperm(size(x,1));
% perm_index2=randperm(size(x,2));
% W = x(perm_index1(1:n),perm_index2(1:N));
% W=(W-repmat(mean(W),size(W,1),1))./repmat(std(W)+1e-8,size(W,1),1);

command = ['/mnt/software/R-3.2.4/bin/Rscript snp_simulation2.R ',num2str(n),' ',num2str(N)];
[status,cmdout] = system(command);
load('/DATA/249/xli/snp2.mat');
W = bsxfun(@rdivide,bsxfun(@minus,snp,mean(snp)),std(snp)+1e-10);

e = normrnd(0, sigmaE, n, 1);
y = W * U + e;

%%
tic
%[u_LHalf, fitinfo] = lasso(W, y);
%parpool()
alpha=1
opts = statset('UseParallel',true);
[u_LHalf, fitinfo] = lasso_coarse(W, y,'CV',10,'Alpha',alpha,'Options',opts);
lassoPlot(u_LHalf,fitinfo,'PlotType','CV');
lambda=fitinfo.Lambda;
fitinfo.IndexMinMSE
fitinfo.Index1SE
toc

%%

% perm_index1=randperm(size(x,1));
% perm_index2=randperm(size(x,2));
% W2 = x(perm_index1(1:n),perm_index2(1:N));
% W2=(W2-repmat(mean(W2),size(W2,1),1))./repmat(std(W2)+1e-8,size(W2,1),1);
% 
% e2 = normrnd(0, sigmaE, n, 1);
% y2 = W2 * U + e2;
% 
% for i=1:length(lambda)-5
%     %[u_LHalf(:,i), cpu_time] = IJT_LHalf(W, y, lambda(i), u_Init, MaxIterNum);
%     %[u_LHalf(:,i), cpu_time] = IJT_L2rds(W, y, lambda(i), u_Init, MaxIterNum);
%     %[u_LHalf(:,i), fitinfo] = lasso(W, y,'lambda',lambda(i));
%     index{i} = find(u_LHalf(:,i));
%     
%     h_cool(i,1) = heritability_cool(y,W(:,index{i}),size(W(:,index{i}),2),size(W(:,index{i}),1));
%     h_cool(i,2) = heritability_cool(y2,W2(:,index{i}),size(W2(:,index{i}),2),size(W2(:,index{i}),1));
%     
%     h_lmm(i,1) = heritability_lmm(y,W(:,index{i}),size(W(:,index{i}),2),size(W(:,index{i}),1));
%     h_lmm(i,2) = heritability_lmm(y2,W2(:,index{i}),size(W2(:,index{i}),2),size(W2(:,index{i}),1));   
% end
% 
% figure
% plot(h_cool(:,1),'-bo');
% hold on 
% plot(h_cool(:,2),'-ro');
% xlabel('lambda')
% ylabel('heritability')
% legend(['trainning set'],['validatoin set'])
% 
% figure
% plot(h_lmm(:,1),'-bo');
% hold on 
% plot(h_lmm(:,2),'-ro');
% xlabel('lambda')
% ylabel('heritability')
% legend(['trainning set'],['validatoin set'])

%%
% A = input('Input a number:');
A=fitinfo.Index1SE;
disp(['You input number is:',num2str(A)]);

% result = index{A};
result = find(u_LHalf(:,A));






