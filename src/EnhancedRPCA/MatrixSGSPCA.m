%% solve the problem 
%              min_X ||D_x(X)||_*+||D_y(X)||_* +||D_z(X)||_* +\lambda||E||_1
%                                   s.t.  Y= X+E
%                          ===============================
%              min_X ||X1||_*+||X2||_* +||X3||_* +\lambda||E||_1
%                            s.t.  Y= X+E
%                                  D_x(X)=X1 
%                                  D_y(X)=X2 
%                                  D_z(X)=X3 
%                          ===============================                       
%         D is difference operator,T is difference tensor,T is known
%  ------------------------------------------------------------------------


function [ output_image] =MatrixSGSPCA(D)
tol     = 1e-6;
maxIter = 40;
rho     = 1.5;
[M,N] = size(D);
max_m = max(M,N);
lambda  = 2/sqrt(max_m);
sizeD   = size(D);
normD   = norm(D,'fro');
% initialize
norm_two = lansvd(D, 1, 'L');
norm_inf = norm( D(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);

mu = 1.25/dual_norm;%1.25/norm_two % this one can be tuned
max_mu = mu * 1e7
%% FFT setting
h               = sizeD(1);
w               = sizeD(2);
%% 
Eny_x   = ( abs(psf2otf([+1; -1], [h,w])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w])) ).^2  ;
determ  =  Eny_x + Eny_y;
%% Initializing optimization variables
X              = randn(M,N);
E              = zeros(M,N);
%M1 =zeros(size(D));  % multiplier for D-X-E
M1 = D / dual_norm;
M2 = M1;%zeros(size(D));  % multiplier for Dx_X-X1
M3 = M2;%zeros(size(D));  % multiplier for Dy_X-X2
% main loop
iter = 0;
tic
while iter<maxIter
    iter          = iter + 1;   
    %% -Updata X1,X2,X3
    [u,s,v] = svd(reshape(diff_x(X,sizeD),[M,N])+M2/mu,'econ');
    X1      = u*softthre(s,1/mu)*v';
    [u,s,v] = svd(reshape(diff_y(X,sizeD),[M,N])+M3/mu,'econ');
    X2      = u*softthre(s,1/mu)*v';
    %% -Updata X
%     x             = X(:);
%     Cha           = D-E;
%     Mu_x          = X1;
%     Mu_y          = X2;
%     Mu_z          = X3;
%     x             = myPCG_sstv(x,Cha,Mu_x,Mu_y,Mu_z,M1,M2,M3,M4,mu,sizeD);
%     X             = reshape(x,[M*N,p]);
    
    diffT_p  = diff_xT(mu*X1-M2,sizeD)+diff_yT(mu*X2-M3,sizeD);
    numer1   = reshape( diffT_p + mu*(D(:)-E(:)) + M1(:), sizeD);
    x        = real( ifftn( fftn(numer1) ./ (mu*determ + mu) ) );
    X        = reshape(x,[M,N]);
%     if mod(iter,10)==0
%         Xten = reshape(X,[M,N,p]);
%         figure;imshow(Xten(:,:,1),[])
%     end
    %% -Update E
    E             = softthre(D-X+M1/mu, lambda/mu);
%     E               = (M1+mu*(D-X))/(2*lambda+mu);% Gaussian noise
    %% stop criterion  
    leq1 = D -X -E;
    leq2 = reshape(diff_x(X,sizeD),[M,N])- X1;
    leq3 = reshape(diff_y(X,sizeD),[M,N])- X2;
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = max(abs(leq2(:)));
    stopC3 = norm(leq3,'fro')/normD;
    disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
            ',Y-X-E=' num2str(stopC1,'%2.3e') ',||DX-X1||=' num2str(stopC2,'%2.3e')...
            ',|DY-X2|' num2str(stopC3,'%2.3e')]);
    if stopC1<tol && stopC2<tol
        break;
    else
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        mu = min(max_mu,mu*rho); 
    end 
%     load('Simu_indian.mat');
%     [mp(iter),sm(iter),er(iter)]=msqia(simu_indian,reshape(X,[M,N,p]));
end
output_image = X;
end