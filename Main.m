clc;
clear all %#ok<CLALL>
% load data
m7g_disease_data = importdata('.\dataset\m7g_disease_data');
wrr5 = importdata('.\dataset\m7GSim');

Wmm = wrr5;
Wdd = importdata('.\dataset\DiseaseSim'); 
share_mm = importdata('.\dataset\share_mm.mat'); 
share_dd = importdata('.\dataset\share_dd.mat'); 
Wrname = importdata('.\dataset\m7gName'); 
Wdname = importdata('.\dataset\DiseasesName.txt');
A_ori= m7g_disease_data;


%% k-fold cross validation
alpha = 0.5;
kfold = 10;
known= find(A_ori);
Indices =crossvalind('Kfold', length(known), kfold);

for h=1:kfold 
test = (Indices == h);
train = ~test;
test_data=known(test,:);
train_data=known(train,:); 

Wmd_train=m7g_disease_data;
Wmd_train(test_data)=0;
Wmd_test = m7g_disease_data;
Wmd_test(train_data)=0;


%% ---cluster-one_1.0----------
[RWmm,RWdd] = Cluster(Wmm,Wdd,share_mm,share_dd,Wrname,Wdname);

%% logistic function
cr = -15;
cd = -15;
dr = log(9999);
dd = log(9999);
LWmm = 1./(1+exp(cr*RWmm+dr)); 
LWdd = 1./(1+exp(cd*RWdd+dd));

%% random walk with restart
Score = RWRnew(Wmd_train,RWmm,RWdd,alpha);
Rt{h} = Score;

%% CV
rt_list{h} = Score(:);
test_list{h} = Wmd_test(:);
m7g_dis = m7g_disease_data;
[TPR_DIP_UP,FPR_DIP_UP,RECALL_DIP_UP,PRECISION_DIP_UP] = ten_CV(Indices,h, m7g_dis,Rt); 
TPR{h} = TPR_DIP_UP;
FPR{h} = FPR_DIP_UP; 
RECALL{h} = RECALL_DIP_UP;
PRECISION{h} = PRECISION_DIP_UP;
end

%% result
TPR_DIP_UP = 0;
FPR_DIP_UP = 0;
RECALL_DIP_UP = 0;
PRECISION_DIP_UP = 0;
for k = 1:kfold
    TPR_DIP_UP=TPR_DIP_UP+TPR{1,k};
    FPR_DIP_UP=FPR_DIP_UP+FPR{1,k};
    RECALL_DIP_UP=RECALL_DIP_UP+RECALL{1,k};
    PRECISION_DIP_UP=PRECISION_DIP_UP+PRECISION{1,k};
end
TPR_DIP_UP = TPR_DIP_UP/kfold;
FPR_DIP_UP = FPR_DIP_UP/kfold;
RECALL_DIP_UP = RECALL_DIP_UP/kfold;
PRECISION_DIP_UP = PRECISION_DIP_UP/kfold;

draw_ROC(FPR_DIP_UP,TPR_DIP_UP,RECALL_DIP_UP,PRECISION_DIP_UP);