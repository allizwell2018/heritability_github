function [h,se] = heritability_lmm(y,W,N,n)
A = W*W'/N ;
[S, D] = eig(A);
lambda = diag(D);

h = ones(1, 1);
h = 0.5*h;
yt = S'*y;

nMaxIter = 1000;
count0  = 0;
count1  = 0;
        
for j = 1:nMaxIter
            
    vec1 = yt.^2.*(lambda-1)./(h*(lambda-1)+1).^2;
    vec2 = yt.^2./(h*(lambda-1)+1);
    vec3 = (lambda-1)./(h*(lambda-1)+1);
    vec4 = yt.^2.*(lambda-1).^2./(h*(lambda-1)+1).^3;
    vec5 = (lambda-1).^2./(h*(lambda-1)+1).^2;          
            
    L1 = ((1/n)*sum(vec1))/((1/n)*sum(vec2))-(1/n)*sum(vec3);
    L2 = ((-2/n)*sum(vec4))/(1/n*sum(vec2))+((1/n)*sum(vec1))^2/((1/n)*sum(vec2))^2+(1/n)*sum(vec5);
    h = h-L1/L2;
            
    if h <= 0
        h = 0.0001;
        count0 = count0+1;
        if count0 == 100
            break
        end
    end
    if h >= 1
        h = 0.9999;
        count1 = count1+1;
        if count1 == 100
            break
        end
    end
           
    if abs(L1) < 1e-6
        break
    end
    if (j == nMaxIter)
        warning('reached %d iterations',nMaxIter)
        break
    end   
end
%% se
Vy = 1;
Vg = h*Vy;
Ve = (1-h)*Vy;
%A = W*W'/m;
V = A*Vg+eye(n)*Ve;
P = inv(V);  
AI = 0.5*[y'*P*A*P*A*P*y y'*P*A*P*P*y; y'*P*P*A*P*y y'*P*P*P*y];
invAI = inv(AI);
VarVg = invAI(1,1);
VarVy = invAI(2,2);
CovVgVe = invAI(1,2);
% VarVy = VarVg+2*CovVgVe+VarVe;
CovVgVy = VarVg+CovVgVe;
Varh = (Vg/Vy)^2*(VarVg/(Vg*Vg)+VarVy/(Vy*Vy)-2*CovVgVy/(Vg*Vy));
se = sqrt(Varh);