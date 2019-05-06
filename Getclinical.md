```
clear 
fileID = fopen('gdcclinicaltrain.txt','r');
slCharacterEncoding('ISO-8859-1')
nNumberCols =14;
format = ['%s ' repmat('%s ', [1 nNumberCols])];
A=textscan(fileID,format,219,'Headerlines',1,'Delimiter','\t');

% Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
% Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
% Current reformed smoker for > 15 years (greater than 15 years) = 3
% Current reformed smoker for ??15 years (less than or equal to 15 years) = 4
% Current reformed smoker, duration not specified = 5
% Smoking History not documented = 7

PatientID=A{1,1};
Tobacco=cellfun(@str2double,A{1,2});
Pind=cellfun(@str2double,A{1,4});
hpid=PatientID(Pind==3);
npid=PatientID(Pind==1);
fclose(fileID);
save luadGclinical.mat PatientID Pind hpid npid

```

