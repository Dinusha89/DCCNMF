clc;
%Parameter 1 -> Dataset name
%Parameter 2 -> Weight Mode : (Binary, HeatKernel, Cosine) values can be used
%Parameter 3 -> Normalization Type : (L1, L2, , MinMax, Std) values can be
%used
MultiRunFile("BBCSport2views.mat", "HeatKernel", "L1");
