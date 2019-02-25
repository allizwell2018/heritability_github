function [h] = heritability_cool2(y,W,d,n)
% sigmaG_cool = -d/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W*W'*y;
% sigmaE_cool = (d+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W*W'*y;

m1 = 1/d*trace(1/n*W'*W);
m2 = 1/d*trace((1/n*W'*W)*(1/n*W'*W))-1/(d*n)*(trace(1/n*W'*W))^2;
sigmaG_cool = -(d*m1^2)/(n*(n+1)*m2)*(y'*y)+m1/(n*(n+1)*m2)*y'*W*W'*y;
sigmaE_cool = (1+(d*m1^2)/((n+1)*m2))*1/n*(y'*y)-m1/(n*(n+1)*m2)*y'*W*W'*y;
h = sigmaG_cool / (sigmaG_cool + sigmaE_cool);
end