%% Low rank matrix decomposition
%% \lambda *N ,s.t. Y = X + N， rank(X)=r
function [clean_image,E] = MSILRMR(Y,lambda,r)
[M,N,p] = size(Y);
D = zeros(M*N,p);
for i=1:p 
  bandp  = Y(:,:,i);
  D(:,i) = bandp(:); 
end
[d,p] = size(D);
%%%%%% 计算D在非零位置的F范数
normD = norm(D, 'fro');
tol = 1e-6;
maxIter = 40;

rho = 1.5;
mu  = 1e-2;
max_mu = 1e6;
[u,s,v] = svd(D,'econ');
X = u(:,1:r)*s(1:r,1:r)*v(:,1:r)';
E = zeros(d,p);
Gam = zeros(d,p);
for iter = 1:maxIter
    Xpre = X;
    %% update the X
    temp    = D-E+Gam/mu;
    [u,s,v] = svd(temp,'econ');
    X       = u(:,1:r)*s(1:r,1:r)*v(:,1:r)';
    %% update the E
%     E       = (mu*(D-X)+Gam)/(mu+2*lambda);
    E       = softthre(D-X+Gam/mu,lambda/mu);
    Gam     = Gam+mu*(D-X-E);
    mu      = min(mu*rho,max_mu);
    %% stop crietial
    err1     = norm(D(:)-X(:)-E(:))/normD;
    err2     = norm(Xpre(:)-X(:))/normD;
    fprintf('the iteration is %d,the error of D-X-E is %f,difference is %f.\n',iter,err1,err2);
    if err1 <=tol  && err2 <=tol
        fprintf('the stop criteria is reached\n');
        break;
    end
end
clean_image = reshape(X,[M,N,p]);
E           = reshape(E,[M,N,p]);
    
