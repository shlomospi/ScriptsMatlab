%plots the bended stracture of 2 rods anclosed by a circle and acted by a force @point B
%we will calculate the function for one rod, rotate it, rotate the other way and connect them to one plot.
%assuming stracture stays the same length, and using the result for n=2 from part one.
%using starting conditions of V(0)=0 and from symmetry V'(d)=0
clc
clear all
r=1;
E=1;
I=1;
d=r*sqrt(2); %rod length

x=0:0.001:d;  %create an x vector for which V will be calculated
V=(x.^3/(6*sqrt(2))-(r.^2/(sqrt(2)))*x)/(E*I); %V as a function of x for named starting conditions
x1=x/sqrt(2)-V/sqrt(2); %rotate x and V 45 degrees
V1=x/sqrt(2)+V/sqrt(2);
x2=2*1.471253493264445-x1; %create a second x vector which will mirror x1 and start from the right side



plot(x1,V1,x2,V1,'b') %plot the structure


