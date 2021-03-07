clear all;
close all;
clc

A = dlmread('SetA.txt');            %Load SetA into Matrix A (2E^5 x 2)
x = A(:,1);                         %Reads first column of A into x
y = A(:,2);                         %Reads second column of A into y   
scatter(x,y)
