function [Rt]=RWR(A_new,RWmm,RWdd,alpha)

normWmm = normFun(RWmm);
normWdd = normFun(RWdd);

R0 = A_new/sum(A_new(:));
Rt=R0;


nRtleft=0;
nRtright=0;
% ---random walk with restart--
l = 4;  % iteration
r = 2;
for t=1:max(l,r)

    if(t<=l)
        nRtleft = alpha * normWmm * Rt + (1-alpha)*R0;
    end
    if(t<=r)
        nRtright = alpha * Rt * normWdd + (1-alpha)*R0;
    end
    Rt =  (ftl*nRtleft + ftr*nRtright)/2);
end
end

