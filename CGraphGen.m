function G = CGraphGen(Outputs,Inputs,th)
%Example:
% Outputs = csvread("InflationOutputs.csv");
% Inputs = csvread("InflationInputs.csv");
% G = CGraphGen(Outputs,Inputs,0.4);
Cm = corr([Outputs Inputs]);
Cm = abs(Cm)>th;
sCm = sum(Cm);
s=[];
t=[];
for k = 1:length(Cm) 
    s = [s k*ones(1,sCm(k))];
    t=[t find(Cm(k,:))];
end
G = digraph(s,t);
plot(G)
end