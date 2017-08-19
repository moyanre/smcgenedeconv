clear all; close all; clc;

% rng('default');
% rng(1000000); %1000


tic;


load('Affy_2cells_Data');
Affy_2cells_Data = double(Affy_2cells_Data);  
pu1 = Affy_2cells_Data(:,31:33); sum_pu1 = sum(pu1'); avg_pu1 = sum_pu1'/3;
pu2 = Affy_2cells_Data(:,1:3); sum_pu2 = sum(pu2'); avg_pu2 = sum_pu2'/3;  
cell_spec_expr_value = [avg_pu1 avg_pu2]; 





        %%% SAMPLES FOR THOROUGH ANALYSIS

        YY1 = Affy_2cells_Data(:,10:12);
        YY2 = Affy_2cells_Data(:,22:24);

        %YY = [YY1 YY2]; % CONTAINS 6 SAMPLES FOR ANALYSIS
        YY = Affy_2cells_Data(:,4:30);
%cell_prop = [0.25 0.25 0.25 0.75 0.75 0.75;...
             %0.75 0.75 0.75 0.25 0.25 0.25];
         
cell_prop = [0.05 0.05 0.05 0.10 0.10 0.10 0.25 0.25 0.25 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.75 0.75 0.75 0.90 0.90 0.90 0.95 0.95 0.95 ;...
             0.95 0.95 0.95 0.90 0.90 0.90 0.75 0.75 0.75 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.25 0.25 0.25 0.10 0.10 0.10 0.05 0.05 0.05];
 
          

ibere = 5000;
M = cell_prop;
X = cell_spec_expr_value(ibere:end,:);
Y = YY(ibere:end,:);

% Y1 = X*M;
% errorr = (Y - Y1);
% err_vec = errorr(:); 
% mean_ML = mean(err_vec);
% varr_global = mean((err_vec - mean_ML).^2); 



I = size(Y,1);
J = size(Y,2);
K = size(M,1);

n_para = I*K + K*J + 1;


% ch1 = 0.0001                ;
% ch2 = 0.0001;
% mid = 0.8;
% 
% E_r1 = 0:ch1:mid-ch1;
% E_r2 = mid:ch2:1;
% E_r = [E_r1 E_r2];
% 
% T = length(E_r);

E_r = 0:0.0002:1;
T = length(E_r);



N = 40; % 10    



%%% PRIOR SPEC
% 1. M

     
mA = 0.5; sdM = 0.1; vA = sdM^2; 
M_bar = mA*ones(K,J);

M_bar_VEC = M_bar(:); 

var_MM = vA*ones(K,J);
var_MM_VEC = var_MM(:);

 
 




% 2. precision
 alpha_prior = 2;
 beta_prior = 100000; %100000




% 3. X

% 
% M_X = [0.4 0.4 0.6 0.6 0.7 0.6;...
%        0.6 0.6 0.4 0.4 0.3 0.4];
% M_X = [M_X];
% 
% for i = 1:I
%     yy = Y(i,:)';
%     DD = M_X';
%     
%     xx = lsqnonneg(DD,yy);
%     X_priors(i,:) = xx';
%     
%     
% end


X_priors = 30*ones(I,K);

var_pri_X = 1000;
vv = 1/var_pri_X;
 
X_priors_tran = X_priors';
X_priors_VECTOR = X_priors_tran(:);
clear X_priors_tran;



wt_all = zeros(T,N);

theta_all = zeros(n_para,N);




%%% KERNEL
   %ksdM = 0.0001;
   %ksdX = 0.2;


for t = 1:T
    t
    T
    
    theta = zeros(n_para,N);
    
    if t == 1
       
        %1. Sample M
          
           for kj = 1:K*J
               
               
            meM = M_bar_VEC(kj);  
            va_M = var_MM_VEC(kj);
            samM = normrnd(meM,sqrt(va_M),[1,N]); 
             
            theta(kj,:) = samM;
             

%                     beM = K*(j-1) + 1;
%                     enM = K*j;
               
           end
         
        
       %2.  Sample lambda i.e precision
   
         theta(K*J + 1,:) = gamrnd(alpha_prior,1/beta_prior,[1 N]);
          
 
           
        %3.  Sample lambda i.e precision
            
            
            for ik = 1:I*K
                
                meX = X_priors_VECTOR(ik);
                sdX = sqrt(var_pri_X);
                
                samX = normrnd(meX,sdX,[1,N]);
                theta(K*J + 1+ ik,:) = samX; 
                
                %theta(K*J + I + ik,:) = samX; 
            end
           
            wt = (1/N)*ones(1,N);
        
       
       
    else  % t = 2,...T
        
        EE = (E_r(t) - E_r(t-1)); 
        
         theta_prev = theta_all;
         
         
         wt_prev = wt_all(t-1,:);
         
         
     
            
            
            
        %%% COMPUTATION OF WEIGHTS
        
         %pdf_lik = zeros(1,N);
         
        for n = 1:N
            
            prev_par = theta_prev(:,n);
            
        
                MM_vec = prev_par(1:K*J);
                MM_mat = vec2mat(MM_vec,K)';
                
                PREC = prev_par(1 + K*J);
                
               
                XX_vec = prev_par(2 + K*J:end);
                %XX_vec = prev_par(1 + K*J + I:end);
                XX_mat = vec2mat(XX_vec,K);
                
                MEAN_Y = XX_mat*MM_mat;
                
                err_sq = (Y - MEAN_Y).^2;
                exp_n = -0.5*EE*(PREC)*sum(sum(err_sq));
                
                pdfY = (PREC^(vpa(I*J*0.5*EE)))*exp(vpa(exp_n));
                
                pdf_lik(n) = pdfY;
           
        end
        
         wt_int = wt_prev.*pdf_lik;
         clear pdf_lik;
        
         wt = wt_int/sum(wt_int);
         wt = double(wt);
         
         
         
         
         
         %%% resampling step 
       
        neff = 1/sum(wt.^2);
        
        if neff<N/10
         
          rr = randsample([1:N],N,true,wt);
    
          theta_prev = theta_prev(:,rr);
          wt = (1/N)*ones(1,N);
        end 
         
        
        
        
        
        
         %%% MUTATION
         
         
         
         for n = 1:N
            
              
             
             thet_lana = theta_prev(:,n);
             
             
             thet_lana_M = thet_lana(1:K*J);
             thet_lana_M_mat = vec2mat(thet_lana_M,K)';
              
             M_M1 = thet_lana_M_mat;
             
             PRECIS = thet_lana(1 + K*J);
             
             thet_lana_X = thet_lana(2+K*J:end);
             thet_lana_X_mat = vec2mat(thet_lana_X,K);
               
             
             
             
             
             %%% move PRECISION
             
             alp_p =  alpha_prior + (I*J*E_r(t))/2;
             beta_p  = beta_prior + E_r(t)*0.5*sum(sum((Y - thet_lana_X_mat*thet_lana_M_mat).^2));
             
             new_prec = gamrnd(alp_p,1/beta_p);
             

             
             E_var = E_r(t)*new_prec;
             
 
              %%% move the M's 
              
              for k = 1:K
                  
                  for j = 1:J
                     
                     meaM = M_bar(k,j);
                     var_M = var_MM(k,j);
                      
                     
                     X_k = thet_lana_X_mat(:,k);
                     Y_j = Y(:,j);
                      
                     u1 = (1/var_M);
                     u2 = E_var*sum(X_k.^2);
                      
                     U_kj = u1 + u2;
                   
                   
                
                   
                   v1 = E_var*sum(Y_j.*X_k);
                   
                   X_del = thet_lana_X_mat;
                   X_del(:,k) = 0;
                   mmj = M_M1(:,j);
                   
                   Bijk =  X_del*mmj;
                   
                   v2 = E_var*sum(X_k.*Bijk);
                   v3 = meaM/var_M;
                   
                   V_kj = v1 - v2 + v3;
                      
                   
                      meNN = V_kj/U_kj;
                      varNN = 1/U_kj;
                      
                      samNN = normrnd(meNN,sqrt(varNN));
                      
                      M_M1(k,j) = samNN;
                      
                  end
                  
              end
             
             
             
              %%% MOVE the X's 
              
              M_here =  M_M1;
              
              
              
              for i = 1:I
                    
                    
                    for k = 1:K
                       
                        mu_ik = X_priors(i,k);
                        M_k = M_here(k,:);
                        Y_i = Y(i,:);
                        
                        
                        
                        
                        
                        %%% Q_ik 
                        
                       
                        
                        u1 = vv;
                        sum_mj_sqr =  sum(M_k.^2);
                        u2 = E_var*sum_mj_sqr;
                      
                        Q_ik = u1 + u2;
                        
                      
                        %%% P_ik 
                        BB1 = E_var*sum(M_k.*Y_i);
                        
                        X_i = thet_lana_X_mat(i,:);
                        X_i(k) = 0; 
                        X_i_minus = X_i;
                        
                        Y_minus = X_i_minus*M_here;
                        BB2 = E_var*sum(M_k.*Y_minus);
                        
                        BB3 = vv*mu_ik;
                        
                        P_ik = BB1 - BB2 + BB3;
                        
                        mea_x_ik = P_ik/Q_ik;
                        var_x_ik = 1/Q_ik;
                        
                        sam_x_ik = normrnd(mea_x_ik,sqrt(var_x_ik));
                        
                        thet_lana_X_mat(i,k) = sam_x_ik;
                       
                    end
                    
              end
              
              
             XXMM = thet_lana_X_mat';
              
             theta(:,n) = [M_here(:); new_prec; XXMM(:)];
              clear XXMM;
             
         end
        
    end
    
    
  wt_all(t,:) = wt;  
  theta_all = theta;
   
  
end


%%% M part

M_part = theta(1:K*J,:);
MB = M_part.*repmat(wt,K*J,1);

mean_M = sum(MB')';

M_est1 = vec2mat(mean_M,K)';
M_est = M_est1./repmat(sum(M_est1),K,1)

%%% prec part
prec_part = theta(1+K*J,:);
ghh = sum(prec_part.*wt);

var_est = 1/ghh

%%% X part
X_part = theta(2+K*J:end,:);
XB = X_part.*repmat(wt,I*K,1);
mean_X = sum(XB')';

X_est = vec2mat(mean_X,K);

cor1 = corr2(X(:,1),X_est(:,1));
cor2 = corr2(X(:,1),X_est(:,2));

if cor1 > cor2
    corr_M = corr2(M(:),M_est(:))
else
    corr_M = -corr2(M(:),M_est(:))
end

time = toc