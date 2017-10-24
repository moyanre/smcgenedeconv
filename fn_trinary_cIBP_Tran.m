function [Bin_Mat_Loni] = fn_trinary_cIBP_Tran(Bin_Mat_Lana,alpha,betas,time_step)



if time_step == 1
    
    K_i = poissrnd(alpha); 
    
        weit = betas/sum(betas);
        first_row = zeros(1,K_i);
        
        for d = 1:K_i
            du = discretesample(weit,1);
            
            if du == 1
               first_row(d) = 0.5;
            else
               first_row(d) = 1;
            end
            
        end

        Bin_Mat_Loni = first_row;

    
    
    
else
   
    non_zero_col = zeros(1,size(Bin_Mat_Lana,2));
    
    for e = 1:size(Bin_Mat_Lana,2)
        col_j = Bin_Mat_Lana(:,e);
        
        if sum(col_j) > 0
           non_zero_col(e) = 1; 
        end
    end
    
    K_plus = sum(non_zero_col);
    
    new_row_part_1 = zeros(1,K_plus);
    
    
    for k = 1:K_plus
        col_k = Bin_Mat_Lana(:,k);
        m_k = length(find(col_k ~= 0));
        
        m_k1 = length(find(col_k == 0.5));
        m_k2 = length(find(col_k == 1));
        
        p1 = (m_k/time_step) * ((betas(1) + m_k1)/(sum(betas) + m_k));
        p2 = (m_k/time_step) * ((betas(2) + m_k2)/(sum(betas) + m_k));
        p0 = 1 - (p1 + p2);
        val = discretesample([p0 p1 p2],1);
        
        if val == 1
            new_row_part_1(k) = 0;  
        elseif val == 2
            new_row_part_1(k) = 0.5;  
        else
            new_row_part_1(k) = 1; 
        end
        
       
        
    end
    
    
    K_i_new = poissrnd(alpha/time_step);
    
            if K_i_new == 0
                
                Bin_Mat_Loni = [Bin_Mat_Lana;new_row_part_1];

            else
                
                
                new_row_part_2 = zeros(1,K_i_new);
                      weit_2 = betas/sum(betas);
                      
                       for v = 1:K_i_new
                           
                           dudu = discretesample(weit_2,1);
                           
                           if dudu == 1
                               new_row_part_2(v) = 0.5;
                           else
                               new_row_part_2(v) = 1;
                           end
                           
                       end
                
                
                new_row = [new_row_part_1 new_row_part_2];


                padding_zeros = zeros(size(Bin_Mat_Lana,1),K_i_new);

                Bin_Mat_Loni = [Bin_Mat_Lana  padding_zeros; new_row];

            end
    
end


end







