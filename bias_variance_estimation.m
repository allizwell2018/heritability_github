function     [bias,variance] = bias_variance_estimation(U,N,n,sigmaE,p,index)

for i=1:100
    for j = 1:N
        W(:, j) = binornd(2, p(j), n, 1);
        %W(:, j) =  W(:,j) ./ sqrt(2*p(j)*(1-p(j)));
        W(:, j) = ( W(:,j) - 2*p(j) ) ./ sqrt(2*p(j)*(1-p(j)));
    end
    
    e = normrnd(0, sigmaE, n, 1);
    %W = normrnd(0, 1, n, N);
    y = W * U + e;
    h_true(i) = ( var(y) - var(e) ) / var(y) ;
    
    
    W2 = W(:,index);
    d2=size(W2,2);
    
    h_estimation(i) = heritability_lmm(y,W2,size(W2,2),size(W2,1));
    
    fprintf('Iteration of %d \n',i)
end

bias = mean(h_true)-mean(h_estimation);
variance = std(h_estimation);

end

