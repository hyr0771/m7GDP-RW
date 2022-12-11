function [ RWmm,RWdd ] = Cluster( Wmm,Wdd,tnWmm,tnWdd,Wrname,Wdname)
dn = size(Wdd,1);
dm = size(Wmm,1);

clc;

crfid = fopen('M7GC.txt','w'); 
fclose(crfid);
diary('M7GC.txt');
diary on
!java -jar "cluster_one-1.0.jar"   "M7GP.txt" -F csv
diary off

cdfid = fopen('DiseasesC.txt','w');
fclose(cdfid);
diary('DiseasesC.txt');
diary on
!java -jar "cluster_one-1.0.jar"  "DiseasesP.txt" -F csv
diary off

%extract cluster results from the M7GC.txt and DiseasesC.txt;

% ClusterM7Gs = cell(100,1);
ClusterM7Gs = cell(100,1);
ClusterM7GsQuality = size(100,1);
crfid = fopen('M7GC.txt','r');
flag = 0;
m7gcluster_n = 0;
while(true)
    tline = fgetl(crfid);
    if(tline==-1)
        break;
    end
    if(flag==1)
        while(true)
            tline = fgetl(crfid);
            if(tline==-1)
                break;
            end
            linestr = regexp(tline,',','split');
            quality = linestr(6);
            pvalue = linestr(7);
            quality = str2double(quality);
            pvalue = str2double(pvalue);

            if(pvalue<0.1)
                m_clusterstr =char(linestr(8));
                m_clusterstr = strtrim(m_clusterstr);
                m_clusterstr(1) = [];
                cstrlen = size(m_clusterstr,2);
                m_clusterstr(cstrlen) = [];
                
                m7gcluster_n = m7gcluster_n+1;
                ClusterM7Gs{m7gcluster_n} = m_clusterstr; 
                ClusterM7GsQuality(m7gcluster_n) = quality;
            end  
        end
    end
    if(flag==1)
        break;
    end
    if(size(strfind(tline,'Detected'),1)~=0)
        flag = 1;
    end
end
fclose(crfid);
ClusterDiseases = cell(100,1);
ClusterDiseasesQuality = size(100,1);
cdfid = fopen('DiseasesC.txt','r');
flag = 0;
diseasecluster_n = 0;
while(true)
    tline = fgetl(cdfid);
    if(tline==-1)
        break;
    end
    if(flag==1)
        while(true)
            tline = fgetl(cdfid);
            if(tline==-1)
                break;
            end
            linestr = regexp(tline,',','split');
            quality = linestr(6);
            pvalue = linestr(7);
            quality = str2double(quality);
            pvalue = str2double(pvalue);
            
            if(pvalue<0.1)
                m_clusterstr =char(linestr(8));
                m_clusterstr = strtrim(m_clusterstr);
                m_clusterstr(1) = [];
                cstrlen = size(m_clusterstr,2);
                m_clusterstr(cstrlen) = [];
                
                diseasecluster_n = diseasecluster_n+1;
                ClusterDiseases{diseasecluster_n} = m_clusterstr;
                ClusterDiseasesQuality(diseasecluster_n) = quality;
            end  
        end
    end
    if(flag==1)
        break;
    end
    if(size(strfind(tline,'Detected'),1)~=0)
        flag = 1;
    end
end
fclose(cdfid);

RWmm = Wmm;
drugPos = zeros(100,1);
%for i = 2:dm
for i = 2:size(tnWmm)
    num = 1;
    cfname = char(Wrname(i));
    for m = 1:m7gcluster_n
        ctline = ClusterM7Gs{m};
        if(size(strfind(ctline,cfname),1)~=0)
            drugPos(num) = m;
            num = num +1;
        end
    end
    if(num>1)
        for j=1:i-1
            Rquality = 0;
            csname = char(Wrname(j));
            tflag = 0;

            for n=1:num-1
                cindex = drugPos(n);
                sctline = ClusterM7Gs{cindex};
                if(size(strfind(sctline,csname),1)~=0)
                    tflag =1;
                    if(Rquality<ClusterM7GsQuality(cindex))
                        Rquality = ClusterM7GsQuality(cindex);
                    end
                end
            end
            if(tflag==1)
                RWmm(i,j) = (1+Rquality) * Wmm(i,j);
                if(RWmm(i,j)>=1)
                    RWmm(i,j) = 0.99;
                end
                if(RWmm(i,j)<Rquality)
                    RWmm(i,j) = Rquality;
                end
                if(RWmm(i,j)>=1)
                    RWmm(i,j) = 0.99;
                end
            end
        end
    end   
end
for i=1:dm
    for j=i+1:dm
        RWmm(i,j) = RWmm(j,i);
    end
end
RWdd = Wdd;
diseasePos = zeros(100,1);
%for i = 2:dn
for i = 2:size(tnWdd)
    num = 1;
    cfname = char(Wdname(i));
    for m = 1:diseasecluster_n
        ctline = ClusterDiseases{m};
        if(size(strfind(ctline,cfname),1)~=0)
            diseasePos(num) = m;
            num = num +1;
        end
    end
    if(num>1)
       for j=1:i-1
            csname = char(Wdname(j));
            tflag = 0;
            Rquality = 0;
            for n=1:num-1
                cindex = diseasePos(n);
                sctline = ClusterDiseases{cindex};
                if(size(strfind(sctline,csname),1)~=0)
                    tflag =1;
                    if(Rquality<ClusterDiseasesQuality(cindex))
                        Rquality = ClusterDiseasesQuality(cindex);
                    end
                end
            end
            if(tflag==1)
                RWdd(i,j) = (1+Rquality)*Wdd(i,j);
                if(RWdd(i,j)>=1)
                    RWdd(i,j) = 0.99;
                end  

                if(RWdd(i,j)<Rquality)
                    RWdd(i,j) = Rquality;
                end
                if(RWdd(i,j)>=1)
                    RWdd(i,j) = 0.99;
                end
            end
        end 
    end
end
for i=1:dn
    for j=i+1:dn
        RWdd(i,j) = RWdd(j,i);
    end
end

for i=1:dm
    RWmm(i,i)=1;
end

for i=1:dn
    RWdd(i,i)=1;
end
end

