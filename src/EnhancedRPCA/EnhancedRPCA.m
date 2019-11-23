%% solve the problem 
%                       min_X
%                       SUM_i(||D_i(X)||_*+\lambda*||D_i(E)||_1),i=x,y,z
%                         s.t. Y = X + E
%                          ===============================
%                          ===============================                       
%         D is difference operator, T is difference tensor,T is known
%  ------------------------------------------------------------------------

function output_image =EnhancedRPCA(oriData3_noise)
tol     = 1e-6;
maxIter = 100;
rho     = 1.5;
% max_mu  = 1e6;
% mu      = 5*1e-2;
% mu1     = 1e-3;
[M,N,p] = size(oriData3_noise);
lambda  = 3/sqrt(M*N);
sizeD   = size(oriData3_noise);
D       = zeros(M*N,p) ;
for i=1:p
    bandp = oriData3_noise(:,:,i);
    D(:,i)= bandp(:);
end
normD   = norm(D,'fro');
% initialize
norm_two = lansvd(D, 1, 'L');
norm_inf = norm( D(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);

mu = 1.25/norm_two % this one can be tuned
mu1 = mu; %*0.1;
max_mu = mu * 1e7
%% FFT setting
h               = sizeD(1);
w               = sizeD(2);
d               = sizeD(3);
%% 
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
determ  =  Eny_x + Eny_y + Eny_z;
%% Initializing optimization variables
X              = randn(M*N,p);
E              = zeros(M*N,p);
X1             = zeros(M*N,p);
X2             = X1;
X3             = X2;
E1             = E;
E2             = E1;
E3             = E2;
%
MD   = D / dual_norm;%zeros(size(D));   % multiplier for D-X-E
MX_1 = zeros(size(D));   % multiplier for \delta(X) - X1
MX_2 = zeros(size(D));   % multiplier for \delta(X) - X2
MX_3 = zeros(size(D));   % multiplier for \delta(X) - X3
ME_1 = zeros(size(D));   % multiplier for \delta(E) - E1
ME_2 = zeros(size(D));   % multiplier for \delta(E) - E2
ME_3 = zeros(size(D));   % multiplier for \delta(E) - E3
% main loop
iter = 0;
tic
while iter<maxIter
    iter          = iter + 1;   
    %% -Update X1 and X2 and X3
    tmp_x         = reshape(diff_x(X,sizeD),[M*N,p])+MX_1/mu;
    [u,s,v]       = svd(tmp_x,'econ');
    X1            = u*softthre(s,1/mu)*v';
    tmp_y         = reshape(diff_y(X,sizeD),[M*N,p])+MX_2/mu;
    [u,s,v]       = svd(tmp_y,'econ');
    X2            = u*softthre(s,1/mu)*v';
    tmp_z         = reshape(diff_z(X,sizeD),[M*N,p])+MX_3/mu;
    [u,s,v]       = svd(tmp_z,'econ');
    X3            = u*softthre(s,1/mu)*v';
    %% -Update E1 and E2 and E3
    E1            = softthre(reshape(diff_x(E,sizeD),[M*N,p])+ME_1/mu1,lambda/mu1);
    E2            = softthre(reshape(diff_y(E,sizeD),[M*N,p])+ME_2/mu1,lambda/mu1);
    E3            = softthre(reshape(diff_z(E,sizeD),[M*N,p])+ME_3/mu1,lambda/mu1);
    %% -Updata X
%     x             = X(:);
%     Cha           = D-E;
%     Mu_x          = X1;
%     Mu_y          = X2;
%     Mu_z          = X3;
%     x             = myPCG_sstv(x,Cha,Mu_x,Mu_y,Mu_z,MD,MX_1,MX_2,MX_3,mu,sizeD);
%     X             = reshape(x,[M*N,p]);
    diffT_p  = diff_xT(mu*X1-MX_1,sizeD)+diff_yT(mu*X2-MX_2,sizeD);
    diffT_p  = diffT_p + diff_zT(mu*X3-MX_3,sizeD);
    numer1   = reshape( diffT_p + mu*(D(:)-E(:)) + MD(:), sizeD);
    x        = real( ifftn( fftn(numer1) ./ (mu*determ + mu) ) );
    X        = reshape(x,[M*N,p]);
    if mod(iter,10)==0
        Xten = reshape(X,[M,N,p]);
        figure;imshow(Xten(:,:,1),[])
    end
    %% -Updata E
%     e             = E(:);
%     Cha           = D-X;
%     Mu_x          = E1;
%     Mu_y          = E2;
%     Mu_z          = E3;
%     e             = myPCG_sstv(e,Cha,Mu_x,Mu_y,Mu_z,MD,ME_1,ME_2,ME_3,mu,sizeD);
%     E             = reshape(e,[M*N,p]);
    diffT_p  = diff_xT(mu1*E1-ME_1,sizeD)+diff_yT(mu1*E2-ME_2,sizeD);
    diffT_p  = diffT_p + diff_zT(mu1*E3-ME_3,sizeD);
    numer1   = reshape( diffT_p + mu*(D(:)-X(:)) + MD(:), sizeD);
    e        = real( ifftn( fftn(numer1) ./ (mu1*determ + mu) ) );
    E        = reshape(e,[M*N,p]);
 
    %% stop criterion  
    leq1 = D -X -E;
    leq2 = reshape(diff_x(X,sizeD),[M*N,p])- X1;
    leq3 = reshape(diff_y(X,sizeD),[M*N,p])- X2;
    leq4 = reshape(diff_z(X,sizeD),[M*N,p])- X3;
    leq5 = reshape(diff_x(E,sizeD),[M*N,p])- E1;
    leq6 = reshape(diff_y(E,sizeD),[M*N,p])- E2;
    leq7 = reshape(diff_z(E,sizeD),[M*N,p])- E3;
    
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = max(abs(leq2(:)));
    stopC4 = norm(leq4,'fro')/normD;
    stopC5 = norm(leq5,'fro')/normD;
    stopC7 = norm(leq7,'fro')/normD;
    disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
            ',Y-X-E=' num2str(stopC1,'%2.3e') ',||DX-UV||=' num2str(stopC2,'%2.3e')...
            ',|DZ-UV|' num2str(stopC4,'%2.3e') ',||D_xE-UV||=' num2str(stopC5,'%2.3e')...
            ',||D_zE-UV||=' num2str(stopC7,'%2.3e')]);
    if stopC1<tol && stopC2<tol
        break;
    else
        MD   = MD   + mu*leq1;
        MX_1 = MX_1 + mu*leq2;
        MX_2 = MX_2 + mu*leq3;
        MX_3 = MX_3 + mu*leq4;
        ME_1 = ME_1 + mu1*leq5;
        ME_2 = ME_2 + mu1*leq6;
        ME_3 = ME_3 + mu1*leq7;
        mu = min(max_mu,mu*rho); 
        mu1 = min(max_mu,mu1*rho);
    end 
%     load('Simu_indian.mat');
%     [mp(iter),sm(iter),er(iter)]=msqia(simu_indian,reshape(X,[M,N,p]));
end
output_image = reshape(X,[M,N,p]);
end