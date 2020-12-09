tic

n = 1000;eps = 1e-1;
x0 = squeeze(ones(1, n) * 3);
X=sym('X', [1, n]);
iter = 0;
dqdrtic_gradient = gradient(dqdrtic(X), X);
A0 = eye(n);
grad0 = double(subs(dqdrtic_gradient, X, x0));
k = 0;
while true
    syms alpha
    dir = (- A0 * dqdrtic_gradient)';
    curdir = double(subs(dir, X, x0));
    iter = iter + 1
    k = k + 1;
    x1 = x0 + alpha * curdir;
    alpha = double( solve(diff(dqdrtic(x1))) );
    x1 = double(x0 + curdir * alpha);
    grad1 = double(subs(dqdrtic_gradient, X, x1));
    sk = double(x1 - x0); yk = double(grad1 - grad0);
    A0 = A0 + sk*(sk')/(sk*yk)-(A0*yk*(yk')*A0)/((yk')*A0*yk);
    x0 = x1;
    grad0 = grad1;
    norm(grad1)
    if norm(grad1) <= eps
        break
    end
    if k == n
        k = 0;
        A0 = eye(n);
        grad0 = double(subs(dqdrtic_gradient, X, x0));
    end
end
x0
dqdrtic(x0)
iter
dir
toc
disp(['运行时间: ',num2str(toc)]);