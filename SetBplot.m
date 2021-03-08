clear all;
close all;
clc

B = dlmread('SetB.txt');            %Load SetA into Matrix A (2E^5 x 2)
xClass1 = B(1:60000,1);             %Reads first column of A, class 1 into x
xClass2 = B(60001:200000,1);        %Reads first column of A, class 2 into x
yClass1 = B(1:60000,2);             %Reads second column of A, class 1 into y
yClass2 = B(60001:200000,2);        %Reads second column of A, class 2 into yx1 = linspace(-10,10,50);
x1 = linspace(-2,4,50);
x2 = linspace(-5,5,50);
x3 = linspace(-2,4,50);
x4 = linspace(-5,5,50);
x2 = sqrt(-.375 .* x1 .* x1 + .5*x2 + 2.886) / .4375;
x4 = -1 * sqrt(-.375 .* x3 .* x3 + .5*x4 + 2.886) / .4375;
scatter(xClass2,yClass2, 'b')
hold on
scatter(xClass1,yClass1, 'g')
hold on
plot(x1,x2,'rd--')
hold on
plot(x3,x4,'rx--')