
% clear
% 
% load height_imagen.mat
% y=height;
% 
% % load pheno_all_imagen.mat
% % subcortical_l=pheno_all(:,[18,19,20,23,26,27,28]);
% % subcortical_r=pheno_all(:,[37,38,39,42,45,46,47]);
% % y=(subcortical_l+subcortical_r)/2;
% 
% % load rhthichness_imagen.mat
% % y=rhthichness;%pheno_all; %y=[zeros(696,1);ones(589,1)];
% 
% 
% gg=load('gene2.txt');
% g=gg(indexr1,:);
% mask=delete_NA(g);
% x=g(:,mask);
% x=(x-repmat(mean(x),size(x,1),1))./repmat(std(x)+1e-8,size(x,1),1);
% 
% y=(y-repmat(mean(y),size(y,1),1))./repmat(std(y)+1e-8,size(y,1),1);

%%
% h=importdata('header2.txt');
% subj=h.textdata(:,1);
% 
% [Nsubj, Ncol, SubjID, Data, PhenoNames] = ParseFile('/DATA/249/xli/973_GWAS/973data/lh.thickness.txt', '\t',0);
% 
% [sub_inter,index1,index2]=intersect(subj,SubjID);

% lambda=[1.228219e-01;1.119107e-01;1.019689e-01;9.291024e-02;8.465635e-02;...
%          7.713572e-02;7.028319e-02;6.403943e-02;5.835034e-02;5.316666e-02;...
%          4.844348e-02;4.413989e-02;4.021863e-02;3.664571e-02;3.339021e-02;...
%          3.042392e-02;2.772114e-02;2.525847e-02;2.301458e-02];

%%
% x=x(index1,:); y=y(index1,:);
% y2=Data(index2,:); y2_old=Data(index2,:); %y2_old2=(y2-repmat(mean(y2),size(y2,1),1));
% y2=(y2-repmat(mean(y2),size(y2,1),1));%./repmat(std(y2),size(y2,1),1);
%parpool(8)
% for i=11:20%size(y2,2)
%     
%  K = x*x'/size(x,2);
%  Cov = ones(size(x,1),1);   
%  heritability(1,i) = heritability_lmm(y2(:,i),x,size(x,2),size(x,1));
%  %heritability(1,i) = heritability_lmm_cov(y2(:,i),x,size(x,2),size(x,1),Cov);
%  %heritability(2,i) = heritability_cool(y2(:,i),x,size(x,2),size(x,1));
%  %[Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat(y2(:,i), Cov, K, 0);
%  [Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat_noCov(y2(:,i), Cov, K, 0);
%  heritability(2,i) =  h2; heritability(3,i) = SE; heritability(4,i) = Pval;
%  
% %opts = statset('UseParallel',true);
% opts = statset('UseParallel',false);
% [u_LHalf, fitinfo] = lasso_coarse(x, y2(:,i),'CV',10,'Alpha',1,'Options',opts);
% lassoPlot(u_LHalf,fitinfo,'PlotType','CV');
% lambda=fitinfo.Lambda;
% fitinfo.IndexMinMSE
% fitinfo.Index1SE
% 
% A = fitinfo.Index1SE;%input('Input a number:');
% %disp(['You input number is:',num2str(A)]);
% 
% index = find(u_LHalf(:,A));
% 
% x2=x(:,index);
% K2 = x2*x2'/size(x2,2);
% heritability(5,i) = heritability_lmm(y2(:,i),x2,size(x2,2),size(x2,1));
% %heritability(5,i) = heritability_lmm_cov(y2(:,i),x2,size(x2,2),size(x2,1),Cov);
% %[Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat(y2(:,i), Cov, K2, 0);
% [Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat_noCov(y2(:,i), Cov, K2, 0);
% heritability(6,i) =  h2; heritability(7,i) = SE; heritability(8,i) = Pval;
% %heritability(4,i) = heritability_cool(y2(:,i),x2,size(x2,2),size(x2,1));
% end
%%

%  K = x*x'/size(x,2);
%  %Cov = ones(size(x,1),1);  
%  [cov,~,~]=svd(K); Cov=cov(:,1:3);
%  %Cov=(Cov-repmat(mean(Cov),size(Cov,1),1));%./repmat(std(y)+1e-8,size(y,1),1);


% randindex=randperm(length(y)); num1=300;
% x1=x(randindex(1:num1),:); y1=y(randindex(1:num1),:);
% x2=x(randindex(num1+1:end),:); y2=y(randindex(num1+1:end),:);

% alpha_list=[0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
alpha_list=[3e-5 1e-4 3e-4 1e-3 3e-3 1];
alpha=alpha_list(4);

 for i=1:size(y,2)

%heritability(1,1) = heritability_lmm(y,x,size(x,2),size(x,1));
[heritability(1,i),heritability(2,i)] = heritability_lmm(y(:,i),x,size(x,2),size(x,1));
%heritability(2,1) = heritability_cool(y(:,i),x,size(x,2),size(x,1));
%[Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat_noCov(y(:,i), Cov, K, 0);
%heritability(3,1) =  h2; heritability(4,1) = SE; heritability(5,1) = Pval;
%[heritability(8,1),heritability(9,1)] = heritability_lmm_cov(y(:,i),x,size(x,2),size(x,1),Cov);%
     
 
opts = statset('UseParallel',true);
[u_LHalf{i}, fitinfo{i}] = lasso_coarse(x1, y1(:,i),'CV',10,'Alpha',alpha,'Options',opts);
%lassoPlot(u_LHalf{i},fitinfo{i},'PlotType','CV');

mse(i)=fitinfo{i}.MSE(fitinfo{i}.IndexMinMSE);
%A=fitinfo{i}.IndexMinMSE;
A=fitinfo{i}.Index1SE;

% A = input('Input a number:');
% disp(['You input number is:',num2str(A)]);

index = find(u_LHalf{i}(:,A));

xx2=x2(:,index);
 
%K2 = xx2*xx2'/size(xx2,2);
%K3=x2*x2'/size(x2,2);
%[cov2,~,~]=svd(K3); Cov2=cov2(:,1:3);
%heritability(1,i+1) = heritability_lmm(y2,xx2,size(xx2,2),size(xx2,1));
[heritability(3,i),heritability(4,i)] = heritability_lmm(y2(:,i),xx2,size(xx2,2),size(xx2,1));
%heritability(2,i) = heritability_cool(y2(:,i),xx2,size(xx2,2),size(xx2,1));%heritability_lmm_cov(y2,xx2,size(xx2,2),size(xx2,1),Cov);%
%[Pval, h2, SE, PermPval, PermFWEcPval, Nsubj, Npheno, Ncov] = MEGHAmat_noCov(y2(:,i), Cov, K2, 0);
%heritability(3,i) =  h2; heritability(4,i) = SE; heritability(5,i) = Pval;
%[heritability(8,i),heritability(9,i)] = heritability_lmm_cov(y2(:,i),xx2,size(xx2,2),size(xx2,1),Cov2);%

 end
 
 %mse
 for j=1:size(y,2)
     A=fitinfo{j}.Index1SE;
     index = find(u_LHalf{j}(:,A));
     heritability(5,j)=length(index)/1e4;
     heritability(6,j)=mse(j);
 end