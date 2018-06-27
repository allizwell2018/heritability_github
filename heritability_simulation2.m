clear
%clc

n_list=[300 1000 3000];
N_list=[100 300 600 1000 3000 6000 10000];

%for nn=1:length(n_list)
%% parameters
n = 1000;  % number of samples
%n = n_list(nn);  % number of samples
N = 10000; % number of snp
d=N;
N_big = 10; % number of big snp effect
N_small = N - N_big;  % number of small snp effect
N_SIS = 9000; % SIS feature reduction
%sigmaY = 1; % total variation
sigmaE = 1; % unexplainable variation
sigmaU_big = sqrt(1/N);%1e-1;
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


for i=1:100
for j = 1:N
    
    Z(:, j) = binornd(2, p(j), n, 1);
    W(:, j) = ( Z(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));   
end
%W = normrnd(0, 1, n, N);
e = normrnd(0, sigmaE, n, 1);
%U = [ normrnd(0, sigmaU_big, N, 1) ];
y = W * U + e;

heritability(i,1) = ( var(y) - var(e) ) / var(y) ;

tic
sigmaG_cool = -d/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W*W'*y;
sigmaE_cool = (d+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W*W'*y;
heritability(i,2) = sigmaG_cool / (sigmaG_cool + sigmaE_cool);
time(i,1)=toc;   

tic
heritability(i,3) = heritability_lmm(y,W,size(W,2),size(W,1));
time(i,2)=toc;

% sigmaG_ols = -d/n/(n-d)*(y'*y) + 1/(n-d)*y'*W/(W'*W) * W'*y;
% sigmaE_ols = 1/(n-d)*(y'*y)-1/(n-d)*y'*W/(W'*W) * W'*y;
% heritability(i,4) = sigmaG_cool / (sigmaG_cool + sigmaE_cool);

fprintf('Iteration of %d \n',i)
end

mean(heritability)
std(heritability)
sum(time)


%end

%%
xnumber=cell(1,3);
xnumber(:,1:3)={'sampling','fix','random'};

xnumber = repmat(xnumber,size(heritability,1),1);
xnumber = xnumber(:);

heritability2=heritability(:);

g=gramm('x',xnumber,'y',heritability2);%,'subset',Lambda~=10);
g.stat_boxplot();
%g.stat_violin('fill','transparent');
g.set_names('x','method ','y','heritability');
g.set_title({'heritability simulation',' ',['number of samples: ',num2str(n)],['number of snp: ',num2str(N)],['sigmaE: ',num2str(sigmaE)]});
g.draw();

