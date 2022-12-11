function [TPR_DIP_UP,FPR_DIP_UP,RECALL_DIP_UP,PRECISION_DIP_UP,TP] = ten_CV(Indices,h, m7g_disease_data,Rt)
test = (Indices == h);
train = ~test;
known= find(m7g_disease_data);
test_data=known(test,:);  
train_data=known(train,:);  

M_recovery_dr_di_cover=Rt{1,h};
M_recovery_dr_di_cover(train_data)=-1;

m7g_disease_data_test=m7g_disease_data;
m7g_disease_data_test(train_data)=-1;

sum_test_associations=sum(m7g_disease_data_test(:)==1);

[m7g_disease_data_rows, m7g_disease_data_cols] = size( m7g_disease_data_test ); 

sort_Rt_data_test=zeros(m7g_disease_data_rows,m7g_disease_data_cols);
index_Rt_data_test=zeros(m7g_disease_data_rows,m7g_disease_data_cols);
for i = 1 : m7g_disease_data_rows   
        [sort_m,idx_m]=sort(M_recovery_dr_di_cover(i,:),'descend'); 
        sort_Rt_data_test(i,:)=sort_m;
        index_Rt_data_test(i,:)=idx_m;       
                           
end



Threshold=(1:1:m7g_disease_data_cols); 

TPR_DIP_UP=zeros(1,length(Threshold));  
FPR_DIP_UP=zeros(1,length(Threshold));  
PRECISION_DIP_UP=zeros(1,length(Threshold));  
RECALL_DIP_UP=zeros(1,length(Threshold)); 

c=zeros(1,m7g_disease_data_rows); 
for a=1:m7g_disease_data_rows
       
        if sum(m7g_disease_data_test(a,:)==1)~=0
            c(a)=a;  
        end
        
end 
sum_test_drug_number=sum(c~=0);


for k=1:length(Threshold) 

   TP=0;
   FN=0;
   FP=0;
   TN=0;
    
   z=0;    

   for a=1:m7g_disease_data_rows
        

        if sum(m7g_disease_data_test(a,:)==1)~=0
           for b=1: Threshold(k)
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==1
                   TP=TP+1;
               end
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==0
                   FP=FP+1;
               end
               
            end   

             for b=Threshold(k)+1: m7g_disease_data_cols
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==1
                   FN=FN+1;
               end
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==0
                   TN=TN+1;
               end

             end
          
             
         end 
      
      
      
     

    if (sum(m7g_disease_data_test(a,:)==1)==0)&&(sum_test_associations-sum_test_drug_number>z)
           
            z=z+1;
          
           for b=1: Threshold(k)
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==1
                   TP=TP+1;
               end
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==0
                   FP=FP+1;
               end

           end   

             for b=Threshold(k)+1: m7g_disease_data_cols
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==1
                   FN=FN+1;
               end
               if m7g_disease_data_test(a,(index_Rt_data_test(a,b)))==0
                   TN=TN+1;
               end

             end
          
             
     end

      
   end


  tpr=TP/(TP+FN);
  fpr=FP/(FP+TN);
  precision=TP/(TP+FP);
  recall=TP/(TP+FN);

  
  TPR_DIP_UP(k)=tpr; 
  FPR_DIP_UP(k)=fpr; 
  PRECISION_DIP_UP(k)=precision;
  RECALL_DIP_UP(k)=recall;
    
   
  
end

end




