function p = compare_two(x,y,nTests)

N = length(x);
all_p = zeros(nTests, 1);
for i=1:nTests
    x = randsample(x,N);
    z = x - y;
    all_p(i) = sum(z>0) / N;
end

p = mean(all_p);
