clc;clear;

%% load data
load lamda_jh_b

N = size(lamda_jh, 1);

lamda_jh = lamda_jh';

Lamda1 = reshape(lamda_jh', 1, N*N);

sum_lamda1 = sum(Lamda1);


Cap_lamda(1, :) = Lamda1 / sum_lamda1;


load lamda_sc_b

N = size(lamda_sc, 1);

lamda_sc = lamda_sc';

Lamda2 = reshape(lamda_sc', 1, N*N);

sum_lamda2 = sum(Lamda2);


Cap_lamda(2, :) = Lamda2 / sum_lamda2;


load lamda_kl_b

N = size(lamda_kl, 1);

lamda_kl= lamda_kl';

Lamda3 = reshape(lamda_kl', 1, N*N);

sum_lamda3 = sum(Lamda3);


Cap_lamda(3, :) = Lamda3 / sum_lamda3;




%% NMF

V = Cap_lamda;

[i, u] = size(V);

K = 10;

W1 = rand(i, K);
B1 = rand(K, u);

maxiter = 1000;
iter = 0;

while 1
    W2 = W1.*((V./(W1*B1))*B1');
%     W2 = W2./(sum(W2,2));
    B2 = B1.*(W2'*(V./(W2*B1)));
    B2 = B2./(sum(B2,2));

    
    iter = iter + 1;
    
    A = W2 * B2;
    
    convcheck = norm(A - V, 'fro') / norm(V, 'fro');
    
    if convcheck < 1e-10 %|| iter == maxiter
        break;
    end
    
    B1 = B2;
    W1 = W2;
    
end
% 
% dim=size(V);                                  
% V=double(V);
% W=rand(dim(1),K);                           
% W=W./(sum(W,2));                   
% 
% B=rand(K,dim(2));
% 
%                                    
% for iter=1:maxiter
%     B=B.*(W'*(V./(W*B)));
%     B=B./(sum(B,2));
%     W=W.*((V./(W*B))*B');
%     
% end

% B2(2:end) = 0;
% B2(1) = 1;

%%
figure
lamda_nmf1 = W2(1, :) * B2;

Lamda_nmf1 = reshape(lamda_nmf1', N, N);

error1 = norm(lamda_nmf1*sum_lamda1 - Lamda1, 'fro');

imagesc(Lamda_nmf1*sum_lamda1); colorbar('FontSize', 15); %colormap(flipud(autumn));
set(gca,'xtick',[], 'ytick', []);
set(gcf,'position',[200,50,800,600])

figure
lamda_nmf2 = W2(2, :) * B2;

Lamda_nmf2 = reshape(lamda_nmf2', N, N);

error2 = norm(lamda_nmf2*sum_lamda2 - Lamda2, 'fro');

imagesc(Lamda_nmf2*sum_lamda2); colorbar('FontSize', 15); %colormap(flipud(autumn));
set(gca,'xtick',[], 'ytick', []);
set(gcf,'position',[200,50,800,600])

figure
lamda_nmf3 = W2(3, :) * B2;

Lamda_nmf3 = reshape(lamda_nmf3', N, N);

error3 = norm(lamda_nmf3*sum_lamda3 - Lamda3, 'fro');

imagesc(Lamda_nmf3*sum_lamda3); colorbar('FontSize', 15); %colormap(flipud(autumn));
set(gca,'xtick',[], 'ytick', []);
set(gcf,'position',[200,50,800,600])

