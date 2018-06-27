%%
function result = lqcv(U,N,n,sigmaE,p)

%lambda = [0.1 1 5 10 15 20 25 30 35 40 45 50 60 70 100];
%lambda=10:5:50;
%lambda = 1e-5:1e-5:1e-4;
MaxIterNum = 1e5;
u_Init = randn(N,1);

for j = 1:N
   
    W(:, j) = binornd(2, p(j), n, 1);
    %W(:, j) =  W(:,j) ./ sqrt(2*p(j)*(1-p(j)));
    W(:, j) = ( W(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
end
e = normrnd(0, sigmaE, n, 1);
y = W * U + e;

for j = 1:N
    
    W2(:, j) = binornd(2, p(j), n, 1);
    %W(:, j) =  W(:,j) ./ sqrt(2*p(j)*(1-p(j)));
    W2(:, j) = ( W2(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
end
e2 = normrnd(0, sigmaE, n, 1);
y2 = W2 * U + e2;

%%

%[u_LHalf, fitinfo] = lasso(W, y);
%parpool()
opts = statset('UseParallel',false);
[u_LHalf, fitinfo] = lasso_coarse(W, y,'CV',10,'Alpha',0.01,'Options',opts);
lassoPlot(u_LHalf,fitinfo,'PlotType','CV');
lambda=fitinfo.Lambda;
fitinfo.IndexMinMSE
fitinfo.Index1SE
for i=1:length(lambda)-5
    %[u_LHalf(:,i), cpu_time] = IJT_LHalf(W, y, lambda(i), u_Init, MaxIterNum);
    %[u_LHalf(:,i), cpu_time] = IJT_L2rds(W, y, lambda(i), u_Init, MaxIterNum);
    %[u_LHalf(:,i), fitinfo] = lasso(W, y,'lambda',lambda(i));
    index{i} = find(u_LHalf(:,i));
    
    h_cool(i,1) = heritability_cool(y,W(:,index{i}),size(W(:,index{i}),2),size(W(:,index{i}),1));
    h_cool(i,2) = heritability_cool(y2,W2(:,index{i}),size(W2(:,index{i}),2),size(W2(:,index{i}),1));
    
    h_lmm(i,1) = heritability_lmm(y,W(:,index{i}),size(W(:,index{i}),2),size(W(:,index{i}),1));
    h_lmm(i,2) = heritability_lmm(y2,W2(:,index{i}),size(W2(:,index{i}),2),size(W2(:,index{i}),1));   
end
lambda2=lambda(1:length(lambda)-5);
figure
semilogx(lambda2,h_cool(:,1),'-bo');
hold on 
semilogx(lambda2,h_cool(:,2),'-ro');
xlabel('Lambda','FontSize',36)
ylabel('Heritability','FontSize',36)
% set(gca,'XTickLabel',lambda)
legend(['Trainning set'],['Validatoin set'])
set(gcf,'Position',[1 1 1539 827])
set(gca,'Position',[0.099 0.252 0.871 0.654])
set(gca,'FontSize',36);


figure
semilogx(lambda2,h_lmm(:,1),'-bo');
hold on 
semilogx(lambda2,h_lmm(:,2),'-ro');
xlabel('Lambda','FontSize',36)
ylabel('Heritability','FontSize',36)
% set(gca,'XTickLabel',lambda)
legend(['Trainning set'],['Validatoin set'])
set(gcf,'Position',[1 1 1539 827])
set(gca,'Position',[0.099 0.252 0.871 0.654])
set(gca,'FontSize',36);



% A = input('Input a number:');
A=fitinfo.Index1SE;
disp(['You input number is:',num2str(A)]);

result = index{A};






