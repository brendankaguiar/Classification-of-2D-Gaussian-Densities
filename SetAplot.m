clear all;
close all;
clc

A = dlmread('SetA.txt');            %Load SetA into Matrix A (2E^5 x 2)
xClass1 = A(1:60000,1);             %Reads first column of A, class 1 into x
xClass2 = A(60001:200000,1);        %Reads first column of A, class 2 into x
yClass1 = A(1:60000,2);             %Reads second column of A, class 1 into y
yClass2 = A(60001:200000,2);        %Reads second column of A, class 2 into y
x1 = linspace(-10,10,50);
x2 = 5.14 - x1;
scatter(xClass1,yClass1, 'g')
hold on
scatter(xClass2,yClass2, 'b')
hold on
plot(x1,x2, 'rd--')

