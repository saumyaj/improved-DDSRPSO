clear all
%mex cec13_func.cpp 
clear all;
close all;
func_num=1;
D=50; % Dimension of the problem
Xmin=-100;
Xmax=100;
pop_size=100; % No of Particles
iter_max=2500; % Max iteration = (10,000 * D)/ps
runs=60; % Number of times to run the algorithm
fhd=str2func('cec13_func'); % Call the test function
fbest = zeros(28,runs);
saved_eval=zeros(28,runs);
fhistory= [];
fileID = fopen('output.txt','w');
a1=1400:-100:100;
a2=-100:-100:-1400;
offset=[a1 a2];
for i=1:28
    func_num=i;
    for j=1:runs
        [gbest,gbestval,fitcount, history,saved]= Diversified_DDSRPSO(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num); 
        fbest(i,j)=gbestval+offset(i);
    end
    s=sort(fbest(i,:));
    f_mean=mean(s(1:30))
    dev=std(fbest(i,:))
    fprintf(fileID,'%d \n',f_mean);
    fprintf(fileID,'%d \n',dev);
end



