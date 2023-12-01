function [H,rearranged_cols]=typical_H(original_H)
%本函数将非典型形式的H转换为典型形式的H。
%对非典型形式H做高斯消去
H_simplified_stair = Gaussian_Elimination(original_H);
%将H化为系统形式
[H_column_permuted, ~,rearranged_cols] = H_systematic(H_simplified_stair, original_H);
%格式转换
H=sparse(logical(H_column_permuted));
