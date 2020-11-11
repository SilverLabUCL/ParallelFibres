



llr = zeros(1,100);

for k = 1:100
    xt = randn + 2.7;
    llr(k) = log(normpdf(xt,3)/normpdf(xt,2.7));
end



    