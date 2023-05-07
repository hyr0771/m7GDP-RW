# m7GDP-RW

We propose a new method named m7GDP-RW to identify m7G-disease associations via random walk with restart. 
This CLTPRW program is implemented in Matlab2018b.

"Main.m" demonstrates the experimental result on the gold standard dataset by m7GDP-RW.

Cluster.m is the function of clusterONE algorithm.

normFun.m is the Laplacian normalization

RWR.m is the random walk with restart algorithm

ten_CV.m is the ten-fold cross-validation

draw_ROC.m used to draw ROC and PR curves

cluster_ One1.0.jar, DiseasesP.txt, M7GP.txt should be placed in the same folder as Main. m

The dataset we used can be found in the "Datasets" folder.
All files of Dataset and Code should be stored in the same folder to run m7GDP-RW.
Contact: 
If you have any questions or suggestions with the code, please let us know. 
Contact Zhihong Wu at wuzh077@126.com or Yiran Huang at hyr@gxu.edu.cn

# About

title={Predicting disease-associated N7‑methylguanosine(m7G) sites via random walk on heterogeneous network }

author={Yiran Huang, Zhihong Wu, Wei Lan, Cheng Zhong*}

year={2022}
