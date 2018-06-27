clear; %close all
%clc
addpath /DATA/249/xli/gramm-master
%% parameters
n = 1000;  % number of samples
N = 10000; % number of snp
d=N;
N_big = 1000; % number of big snp effect
N_small = N - N_big;  % number of small snp effect
N_SIS = 9000; % SIS feature reduction
%sigmaY = 1; % total variation
sigmaE = 1; % unexplainable variation
sigmaU_big = sqrt(1/N_big);%1e-1;
sigmaU_small = 0;
Lambda = 50;
%heritability_gcta = (N_big*sigmaU_big+N_small*sigmaU_small)/(N_big*sigmaU_big+N_small*sigmaU_small+sigmaE);
heritability_gcta = (N*sigmaU_big^2)/(N*sigmaU_big^2+sigmaE^2);
%% 

p = unifrnd(0.1, 0.5 , 1, N);
W = zeros(n, N);
%U = [ normrnd(0, sigmaU_big, N, 1) ];
U = [ normrnd(0, sigmaU_big, N_big, 1) ; normrnd(0, sigmaU_small, N_small, 1) ];
%U = [ 0.1+unifrnd(0.1, 0.5, N_big, 1) ; normrnd(0, sigmaU_small, N_small, 1) ];
%U = [ 1.1;0.9;-1;0.6;0.3;-1.4;1.3;1.8;-2.4;-1.3; normrnd(0, sigmaU_small, N_small, 1) ];

% perm = randperm(N);
% U = U(perm);

%%
index = [1:N_big];
index_lq = lqcv(U,N,n,sigmaE,p);
%index_lq = siscv(U,N,n,sigmaE,p);

for i=1:100
for j = 1:N
    W(:, j) = binornd(2, p(j), n, 1);
    %W(:, j) =  W(:,j) ./ sqrt(2*p(j)*(1-p(j)));
    W(:, j) = ( W(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
end
% for j = N_big:N
%     p(j) = unifrnd(0.1, 0.5);
%     W(:, j) = binornd(2, p(j), n, 1);
%     W(:, j) =  W(:,j)  ./ sqrt(2*p(j)*(1-p(j)));
%     W(:, j) = ( W(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
% end
e = normrnd(0, sigmaE, n, 1);
%W = normrnd(0, 1, n, N);
y = W * U + e;
heritability(i,1) = ( var(y) - var(e) ) / var(y) ;
%  figure
%  [f1,xi1]=ksdensity(y);
%  plot(xi1,f1,'b')
%  hold on
%  [f2,xi2]=ksdensity(W*U);
%  plot(xi2,f2,'r')

W2 = W(:,index);
d2=size(W2,2);

W_lq = W(:,index_lq);

tic
sigmaG_cool = -d/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W*W'*y;
sigmaE_cool = (d+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W*W'*y;
heritability(i,2) = sigmaG_cool / (sigmaG_cool + sigmaE_cool);
time(i,1)=toc;   

tic
heritability(i,3) = heritability_lmm(y,W,size(W,2),size(W,1));
time(i,2)=toc;

sigmaG_cool2 = -d2/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W2*W2'*y;
sigmaE_cool2 = (d2+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W2*W2'*y;
heritability(i,4) = sigmaG_cool2 / (sigmaG_cool2 + sigmaE_cool2);

heritability(i,5) = heritability_lmm(y,W2,size(W2,2),size(W2,1));

heritability(i,6) = heritability_cool(y,W_lq,size(W_lq,2),size(W_lq,1));
heritability(i,7) = heritability_lmm(y,W_lq,size(W_lq,2),size(W_lq,1));

fprintf('Iteration of %d \n',i)
end

figure
% boxplot(heritability,{'sampling','no_sparse_cool','no_sparse_lmm','true_sparse_cool','true_sparse_lmm','lq_sparse_cool','lq_sparse_lmm'})
% ylabel('heritability')

xnumber=cell(1,7);
xnumber(:,1:7)={'sampling','fix','random','fix_sparse_true','random_sparse_true','fix_sparse_lq','random_sparse_lq'};

xnumber = repmat(xnumber,size(heritability,1),1);
xnumber = xnumber(:);

heritability2=heritability(:);

g=gramm('x',xnumber,'y',heritability2);%,'subset',Lambda~=10);
g.stat_boxplot();
%g.stat_violin('fill','transparent');
g.set_names('x',' ','y',' ');
g.set_title({'heritability simulation',' ',['number of samples: ',num2str(n)],['number of snp: ',num2str(N)],['number of sparse snp: ',num2str(N_big)]});
g.draw();


% ylabel('Estimated Heritability','FontName','Helvetica','FontSize',24)
% title({['(a) Number of Causal SNPs: s=',num2str(N_big)]},'FontSize',24);
% axis([ 0 8 0 1])
% box on
% set(gca,'XTickLabel',{'fixed','fixed_SpaR','fixed_oracle','rand.','rand._SpaR','rand._oracle','approx.'})
% set(gca,'FontSize',20);
% set(gcf,'Position',[1 1 1410 591])
% % set(gca,'Position',[0.08 0.09 0.862 0.676])
% set(gca,'Position',[0.099 0.166 0.871 0.724])
% xlabel_rotation

mean(heritability)
std(heritability)
sum(time)
