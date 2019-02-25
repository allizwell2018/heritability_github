%%
function result = lqcv_bias_variance(U,N,n,sigmaE,p)

% dbstop at 64
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


%%

%[u_LHalf, fitinfo] = lasso(W, y);
%parpool()
opts = statset('UseParallel',true);
[u_LHalf, fitinfo] = lasso_coarse(W, y,'CV',10,'Alpha',1,'Options',opts);
lassoPlot(u_LHalf,fitinfo,'PlotType','CV');
lambda=fitinfo.Lambda;
fitinfo.IndexMinMSE
fitinfo.Index1SE
length(lambda)
for i=1:length(lambda)-5
    
    index{i} = find(u_LHalf(:,i));

    [bias(i),variance(i)] = bias_variance_estimation(U,N,n,sigmaE,p,index{i});
    
end
lambda2=lambda(1:length(lambda)-5);
figure

% semilogx(lambda2,bias,'-bo');
% hold on 
% semilogx(lambda2,variance,'-ro');
% xlabel('Lambda','FontSize',36)
% ylabel('Heritability','FontSize',36)
% % set(gca,'XTickLabel',lambda)
% legend(['Bias'],['Variance'])

[hAx,hLine1,hLine2]=plotyy(lambda2,bias,lambda2,variance,'semilogx','semilogx');
% [hAx,hLine1,hLine2]=plotyy(lambda2(end:-1:1),bias(end:-1:1),lambda2(end:-1:1),variance(end:-1:1),'semilogx','semilogx');
xlabel('Lambda','FontSize',36)
ylabel(hAx(1),'Bias','FontSize',36)
ylabel(hAx(2),'Variance','FontSize',36)
set(hLine1,'LineStyle','--','LineWidth', 1)
set(hLine2,'LineStyle',':','LineWidth', 1)

lambda(fitinfo.IndexMinMSE)
lambda(fitinfo.Index1SE)
line([lambda(fitinfo.IndexMinMSE) lambda(fitinfo.IndexMinMSE)] ,[0 0.4],'LineStyle',':','LineWidth', 1)
line([lambda(fitinfo.Index1SE) lambda(fitinfo.Index1SE)] ,[0 0.4],'LineStyle',':','LineWidth', 1)


% set(gca,'Xdir','reverse') set(gca,'Xdir','normal')

% set(gcf,'Position',[1 1 1539 827])
% set(gca,'Position',[0.099 0.252 0.871 0.654])
% set(gca,'FontSize',36);


A = input('Input a number:');
% A=fitinfo.Index1SE;
disp(['You input number is:',num2str(A)]);

result = index{A};






