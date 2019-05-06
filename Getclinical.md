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

------

<!--Get mutation of LUAD & LUSC from GDC-->

```
clear all
fileID = fopen('gdcmutect.maf','r');
slCharacterEncoding('ISO-8859-1')
nNumberCols =120;
format = ['%s ' repmat('%s ', [1 nNumberCols])];
A=textscan(fileID,format,208186,'Headerlines', 6,'Delimiter','\t');
%keyboard
HugoSymbol=A{1,1};
EntrezGene=A{1,2};
chrom=A{1,5};
startpo=A{1,6};
endpo=A{1,7};
VariantClass=A{1,9};
VariantType=A{1,10};
Ref=A{1,11};
Tum1=A{1,12};
Tum2=A{1,13};
pidtemp=char(A{1,16});
pidtemp2=cellstr(pidtemp(:,14:16));
is01A=ismember(pidtemp2,{'01A'})~=0;
HugoSymbol=HugoSymbol(is01A);
VariantClass=VariantClass(is01A);
VariantType=VariantType(is01A);
Tum1=Tum1(is01A);
Tum2=Tum2(is01A);
pidt2=cellstr(pidtemp(:,1:12));
pidt=pidt2(is01A);
fclose(fileID);

save luadGmu.mat HugoSymbol EntrezGene VariantClass VariantType Tum1 Tum2 pidt
```

