tic
n = 10;eps = 1e-1;
x0 = squeeze(ones(1, n) * 3);
X=sym('X', [1, n]);
iter = 0;
syms X
for iter = 1:100000
    d = mod(iter, n) + 1;
    iter
    x1 = [x0(1:d-1), X, x0(d+1:n)];
    dqdrtic_diff = diff(dqdrtic(x1), X);
    x1(d) = double(solve(dqdrtic_diff));
    if norm(x1-x0) <= eps
        break
    end
    x0 = x1;
end
x0
dqdrtic(x0)
iter

toc
disp(['运行时间: ',num2str(toc)]);