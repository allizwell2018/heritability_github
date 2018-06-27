function [h] = heritability_cool(y,W,d,n)
sigmaG_cool = -d/(n*(n+1))*(y'*y)+1/(n*(n+1))*y'*W*W'*y;
sigmaE_cool = (d+n+1)/(n*(n+1))*(y'*y)-1/(n*(n+1))*y'*W*W'*y;
h = sigmaG_cool / (sigmaG_cool + sigmaE_cool);
end