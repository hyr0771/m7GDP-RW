function [Rt]=RWRnew(A_new,RWrr,RWdd,alpha)

normWrr = normFun(RWrr);
normWdd = normFun(RWdd);

R0 = A_new/sum(A_new(:));
Rt=R0;


nRtleft=0;
nRtright=0;
% ---random walk with restart--------------------------------------------
l = 4;   % 随机游走的最大迭代次数
r = 2;
count = 0;
residue = 1;
threeshold = 0.00001;
while (residue >= threeshold)
    R_old = Rt;
    for t=1:max(l,r)

        if(t<=l)
            nRtleft = alpha * normWrr * Rt + (1-alpha)*R0;
            ftl = t;
        end
        if(t<=r)
            nRtright = alpha * Rt * normWdd + (1-alpha)*R0;
            ftr = t;
        end
        Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
    end
    residue = sqrt(sum((Rt - R_old).^2));
    count = count + 1;
end
end

