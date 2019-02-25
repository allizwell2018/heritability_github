clear
%clc
addpath /DATA/249/xli/gramm-master
n_list=[300  1000  3000];
N_list=[100 300 1000 3000 10000]; N_list=[1000 3000 10000];
h_list=[0.25 1 4];
list=N_list;
for nn=1:length(list) 
%% parameters
n = 1000;  % number of samples
% n = list(nn);  % number of samples
% N = 10000; % number of snp
N = list(nn);
d=N;
N_big = 10; % number of big snp effect
N_small = N - N_big;  % number of small snp effect
N_SIS = 9000; % SIS feature reduction
%sigmaY = 1; % total variation
sigmaE = 1; % unexplainable variation
sigmaU_big = sqrt(1/N);%1e-1;
% sigmaU_big = sqrt(h_list(nn)/N);%1e-1;
sigmaU_small = 0;
Lambda = 50;
%heritability_gcta = (N_big*sigmaU_big+N_small*sigmaU_small)/(N_big*sigmaU_big+N_small*sigmaU_small+sigmaE);
heritability_gcta = (N*sigmaU_big^2)/(N*sigmaU_big^2+sigmaE^2);

%%

p = unifrnd(0.1, 0.5 , 1, N);
Z = zeros(n, N);
W = zeros(n, N);
U = [ normrnd(0, sigmaU_big, N, 1) ];
%U = [ unifrnd(0, 0.1, N, 1)] ;
%U = [ normrnd(0, sigmaU_big, N_big, 1) ; normrnd(0, sigmaU_small, N_small, 1) ];
%U = [ 1.1;0.9;-1;0.6;0.3;-1.4;1.3;1.8;-2.4;-1.3;0.7 ; normrnd(0, sigmaU_small, N_small, 1) ];

% perm = randperm(N);
% U = U(perm);


for i=1:50
i   
% for j = 1:N   
%     Z(:, j) = binornd(2, p(j), n, 1);
%     W(:, j) = ( Z(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
% end

% command = '/mnt/software/R-3.2.4/bin/Rscript snp_simulation.R';
command = ['/mnt/software/R-3.2.4/bin/Rscript snp_simulation.R ',num2str(n),' ',num2str(N)];
[status,cmdout] = system(command);
load('/DATA/249/xli/snp.mat');
W = bsxfun(@rdivide,bsxfun(@minus,snp,mean(snp)),std(snp)+1e-10);

%W = normrnd(0, 1, n, N);
e = normrnd(0, sigmaE, n, 1);
U = [ normrnd(0, sigmaU_big, N, 1) ];
y = W * U + e;

heritability(i,nn,1) = ( var(y) - var(e) ) / var(y) ;

tic
% sigmaG_cool = -d/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W*W'*y;
% sigmaE_cool = (d+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W*W'*y;
% heritability(i,nn,2) = sigmaG_cool / (sigmaG_cool + sigmaE_cool);

heritability(i,nn,2) = heritability_cool2(y,W,size(W,2),size(W,1));
time(i,1)=toc;   

tic
heritability(i,nn,3) = heritability_lmm(y,W,size(W,2),size(W,1));
time(i,2)=toc;

% sigmaG_ols = -d/n/(n-d)*(y'*y) + 1/(n-d)*y'*W/(W'*W) * W'*y;
% sigmaE_ols = 1/(n-d)*(y'*y)-1/(n-d)*y'*W/(W'*W) * W'*y;
% heritability(i,4) = sigmaG_cool / (sigmaG_cool + sigmaE_cool);

%fprintf('Iteration of %d \n',i)
end

% mean(heritability)
% std(heritability)
% sum(time)

end

xnumber = repmat(list,size(heritability,1),1,size(heritability,3));
xnumber = xnumber(:);
xnumber = num2str(xnumber,'%3.0f');
xnumber =  mat2cell(xnumber,ones(size(xnumber,1),1),size(xnumber,2));

Lambda=cell(1,1,3);
Lambda(:,:,1:3)={'sampling','fix','random'};

Lambda = repmat(Lambda,size(heritability,1),size(heritability,2),1);
Lambda = Lambda(:);

heritability2=heritability(:);

% boxplot(error2,{sparsity,Lambda},'colorgroup',Lambda)
g=gramm('x',xnumber,'y',heritability2,'color',Lambda);%,'subset',Lambda~=10);
g.stat_boxplot();
%g.stat_violin('fill','transparent');
g.set_names('x','number of sample','y','heritability','color','method');
g.set_title({'heritability simulation',' ','number of snp: 10000'});
g.draw();

if 0
box on
axis([0 4 0 1])
ylabel('Estimated Heritability','FontName','Helvetica','FontSize',28)
% set(gca,'XTickLabel',{'fixed','rand.','approx.'})
% set(gca,'XTickLabel',{'4','1','0.25'})
set(gca,'FontSize',24);
set(gcf,'Position',[1 1 1539 827])
set(gca,'Position',[0.099 0.252 0.871 0.654])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

print( gcf,'-r300','-depsc','figure/NumberOfSamples');
print( gcf,'-r300','-dpng','figure/NumberOfSamples');
saveas(gcf,'figure/NumberOfSamples.fig');
end
print( gcf,'-r300','-depsc','figure/NoiseLevel');
print( gcf,'-r300','-dpng','figure/NoiseLevel');
saveas(gcf,'figure/NoiseLevel.fig');

print( gcf,'-r300','-depsc','figure/NumberOfCausal');
print( gcf,'-r300','-dpng','figure/NumberOfCausal');
saveas(gcf,'figure/NumberOfCausal.fig');

print( gcf,'-r300','-depsc','figure/NumberOfSNPs');
print( gcf,'-r300','-dpng','figure/NumberOfSNPs');
saveas(gcf,'figure/NumberOfSNPs.fig');

