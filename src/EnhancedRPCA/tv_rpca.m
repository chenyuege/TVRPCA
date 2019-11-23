%% solve the problem 
%              min_X ||D_x(X)||_*+||D_y(X)||_* +||D_z(X)||_* +\lambda||E||_1
%                                   s.t.  Y= X+E
%              if there is no smoothless property in one's dimension
%                 there is no regularization term
%                          ===============================
%              min_X ||X1||_*+||X2||_* +||X3||_* +\lambda||E||_1
%                            s.t.  Y= X+E, X is matrix
%                                  D_x(X)=X1 
%                                  D_y(X)=X2 
%                                  D_z(X)=X3 
%                          ===============================                       
%         D is difference operator,T is difference tensor,T is known
%  ------------------------------------------------------------------------

function [A_hat,E_hat] =tv_rpca(NoiseMatrix,tensorsize,tv_flag)

[M,p]=size(NoiseMatrix);
if (tensorsize(1)*tensorsize(2) ~=M) &&(tensorsize(3) ~=p)
    disp('the size of tensor is wrong, please adjust the size parameter')
end
lambda   = sum(tv_flag)/sqrt(max(M,p));
disp(['   the tv regularization term number is :=   ',num2str(sum(tv_flag))]);
if lambda ==0
    lambda = 1/sqrt(max(M,p));
    tol    = 1e-6;
    maxIter=100;
    [A_hat,~,~] = inexact_alm_rpca(NoiseMatrix, lambda, tol, maxIter);
    E_hat = NoiseMatrix - A_hat;
else
    tol     = 1e-6;
    maxIter = 100;
    rho     = 1.5;
    normD    = norm(NoiseMatrix,'fro');
    % initialize
    norm_two  = lansvd(NoiseMatrix, 1, 'L');
    norm_inf  = norm( NoiseMatrix(:), inf) / lambda;
    dual_norm = max(norm_two, norm_inf);
    mu        = 1.25/dual_norm;%1.25/norm_two % this one can be tuned
    max_mu    = mu * 1e7;
    disp(['    mu = ',num2str(mu),'    max_mu = ',num2str(max_mu),' ... ']);
    %% FFT setting
    h               = tensorsize(1);
    w               = tensorsize(2);
    d               = tensorsize(3);
    %% 
    Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
    Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
    Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
    Eny_z   =  permute(Eny_z, [3, 1 2]);
    determ  =  tv_flag(1)*Eny_x + tv_flag(2)*Eny_y + tv_flag(3)*Eny_z;
    %% Initializing optimization variables
    A_hat              = randn(M,p);
    E_hat              = zeros(M,p);
    %M1 =zeros(size(D));  % multiplier for D-X-E
    M1 = NoiseMatrix / dual_norm;
    M2 = M1;%zeros(size(D));  % multiplier for Dx_X-X1
    M3 = M2;%zeros(size(D));  % multiplier for Dy_X-X2
    M4 = M3;%zeros(size(D));  % multiplier for Dz_X-X3
    % main loop
    iter = 0;
    tic
    while iter<maxIter
        iter          = iter + 1;   
        %% -Updata X1,X2,X3
        if tv_flag(1) ==1
            [u,s,v] = svd(reshape(diff_x(A_hat,tensorsize),[M,p])+M2/mu,'econ');
            X1      = u*softthre(s,1/mu)*v';
        else
            X1      = zeros(M,p);
        end
        if tv_flag(2) ==1
            [u,s,v] = svd(reshape(diff_y(A_hat,tensorsize),[M,p])+M3/mu,'econ');
            X2      = u*softthre(s,1/mu)*v';
        else
            X2      = zeros(M,p);
        end
        if tv_flag(3) ==1
            [u,s,v] = svd(reshape(diff_z(A_hat,tensorsize),[M,p])+M4/mu,'econ');
            X3      = u*softthre(s,1/mu)*v';
        else
            X3      = zeros(M,p);
        end
        %% -Updata X
        diffT_p  = tv_flag(1)*diff_xT(mu*X1-M2,tensorsize);
        diffT_p  = diffT_p + tv_flag(2)*diff_yT(mu*X2-M3,tensorsize);
        diffT_p  = diffT_p + tv_flag(3)*diff_zT(mu*X3-M4,tensorsize);
        numer1   = reshape( diffT_p + mu*(NoiseMatrix(:)-E_hat(:)) + M1(:), tensorsize);
        x        = real( ifftn( fftn(numer1) ./ (mu*determ + mu) ) );
        A_hat        = reshape(x,[M,p]);
        %% -Update E
        E_hat             = softthre(NoiseMatrix-A_hat+M1/mu, lambda/mu);
    %     E               = (M1+mu*(D-X))/(2*lambda+mu);% Gaussian noise
        %% stop criterion  
        leq1 = NoiseMatrix -A_hat -E_hat;
        if tv_flag(1) ==1
            leq2 = reshape(diff_x(A_hat,tensorsize),[M,p])- X1;
        else
            leq2 = zeros(M,p);
        end
        if tv_flag(2) ==1
            leq3 = reshape(diff_y(A_hat,tensorsize),[M,p])- X2;
        else
            leq3 = zeros(M,p);
        end
        if tv_flag(3) ==1
            leq4 = reshape(diff_z(A_hat,tensorsize),[M,p])- X3;
        else
            leq4 = zeros(M,p);
        end
        stopC1 = norm(leq1,'fro')/normD;
        stopC2 = max([max(abs(leq2(:))),max(abs(leq3(:))),max(abs(leq4(:)))]);
        stopC3 = max([norm(leq2,'fro'),norm(leq3,'fro'),norm(leq4,'fro')])/normD;
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
                ',Y-X-E=' num2str(stopC1,'%2.4e') ',max norm_1 error =' num2str(stopC2,'%2.4e')...
                ',max norm_fro error =' num2str(stopC3,'%2.4e')]);
        if stopC1<tol && stopC2<tol
            break;
        else
            M1 = M1 + mu*leq1;
            M2 = M2 + mu*leq2;
            M3 = M3 + mu*leq3;
            M4 = M4 + mu*leq4;
            mu = min(max_mu,mu*rho); 
        end 
    %     load('Simu_indian.mat');
    %     [mp(iter),sm(iter),er(iter)]=msqia(simu_indian,reshape(X,[M,N,p]));
    end
end
end