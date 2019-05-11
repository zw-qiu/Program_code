> #### Get clinical information of  tumors

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

> #### Get mutation of LUAD & LUSC from GDC

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

------

> Get TMB /TTR /INDEL-Scores of each gene

```
clear
load luadGmu.mat HugoSymbol VariantType Tum1 Tum2 pidt %Ref EntrezGene chrom startpo endpo VariantClass 
load luadGclinical.mat  PatientID Pind

ish=find(Pind==3);
isn=find(Pind==1);
pidh=PatientID(ish);
pidn=PatientID(isn);

caseid=unique(pidt,'stable');%所有病人名称
geneall=unique(HugoSymbol,'stable');
base={'C';'G';'T';'A'};
[geneid,~,~]=unique(HugoSymbol,'stable');%所有基因名称
[q1,~]=ismember(pidn,pidt);
[q2,~]=ismember(pidh,pidt);
pidn(q1==0)=[];%不吸烟患者的病人名称 
pidh(q2==0)=[];%吸烟患者的病人名称 
pid=[pidh;pidn];
Pind=[3*ones(length(pidh),1);ones(length(pidn),1)];
save patient.mat pidh pidn pid Pind

%% heavy&never ecah gene allmucount(including SNV DEL INS)
geneallmu=zeros(length(Pind),length(geneall));%heavy
for i=1:length(pid)
    %吸烟患者突变信息
    isit=find(ismember(pidt,pid(i))~=0);
    geneaut=HugoSymbol(isit);
    varaut=VariantType(isit);
    tum1aut=Tum1(isit);
    tum2aut=Tum2(isit);
    [unigene,ia,~] = unique(geneaut,'stable');
    for m=1:length(ia)
    genemu(m)=length(find(ismember(geneaut,unigene(m))~=0));%每个病人基因突变的个数
    end
    [~,loc]=ismember(unigene,geneall);
    geneallmu(i,loc)=genemu;
    clear unigene ia genemu loc
end

%% heavy and never each gene TTratio & Del Ins vscore
genealltt=zeros(length(Pind),length(geneall));
genedel=zeros(length(Pind),length(geneall));
geneins=zeros(length(Pind),length(geneall));
for i=1:length(pid)
    %吸烟患者突变信息
    isit=find(ismember(pidt,pid(i))~=0);
    geneaut=HugoSymbol(isit);
    varaut=VariantType(isit);
    tum1aut=Tum1(isit);
    tum2aut=Tum2(isit);
    issnp=find(ismember(varaut,{'SNP'})~=0);
    genenow1=geneaut(issnp);
    tum1now1=tum1aut(issnp);
    tum2now1=tum2aut(issnp);
    [unigene2,ib,~] = unique(genenow1,'stable');%C-unique gene; ia-unique index；ic-genenow1 index
    for n=1:length(ib)%each gene
       isg1=find(ismember(genenow1,unigene2(n))~=0);
       tum1now2=tum1now1(isg1);
       tum2now2=tum2now1(isg1);
       
       count=zeros(4,4);
    for j=1:4
        isA=ismember(tum1now2,base{j});
        tum2now11=tum2now2;
        tum2temp=tum2now11(isA);
        tt=1:4;
        tt(j)=[];
        for k=tt
            if ~isempty(tum2temp)
                isit2=ismember(tum2temp,base{k});
                count(j,k)=length(find(isit2~=0));
                clear isit2
            else
                count(j,k)=0;
            end
        end
        clear tum2temp
    end
    %transversin
    ca(n)=count(1,4);
    gt(n)=count(2,3);
    %transition
    ct(n)=count(1,3);
    ga(n)=count(2,4);
    %ta(l)=count(3,4);
    %tg(l)=count(3,2);
    %tc(l)=count(3,1);
    %gc(l)=count(2,1);
    %cg(l)=count(1,2); 
    %at(l)=count(4,3);
    %ac(l)=count(4,1);
    %ag(l)=count(4,2);
    transversion(n)=ca(n)+gt(n);
    transition(n)=ct(n)+ga(n);
    if transition(n)~=0
        genettratio(n)=transversion(n)/transition(n);
    else
        genettratio(n)=0;
    end
    end
    [~,loc2]=ismember(unigene2,geneall);
    genealltt(i,loc2)=genettratio;
    clear unigene2 genettratio loc2 genettratio
     %突变为DEL vscore
    isdel=find(ismember(varaut,{'DEL'})~=0);
    genenow3=geneaut(isdel);
    tum1now3=tum1aut(isdel);
    tum2now3=tum2aut(isdel);
    [unigene3,ic,~] = unique(genenow3,'stable');
    for n3=1:length(ic)
        isg2=find(ismember(genenow3,unigene3(n3))~=0);
        tum1now22=tum1now3{isg2};
        countc=length(find(tum1now22=='C'));
        countg=length(find(tum1now22=='G'));
        countt=length(find(tum1now22=='T'));
        counta=length(find(tum1now22=='A'));
        xd(n3)=(countc+countt)-(counta+countg);
        yd(n3)=(countg+countt)-(counta+countc);
        zd(n3)=(countg+countc)-(counta+countt);         
    end
    if n3~=0
        delscore=(xd.^2+yd.^2+zd.^2+1)/4;
        [~,loc3]=ismember(unigene3,geneall);
    else
        delscore=zeros(1,length(geneall));
        loc3=1:length(geneall);
    end
    genedel(i,loc3)=delscore;
    clear unigene3 delscore loc3 xd yd zd
    %突变为INS vscore
    isins=find(ismember(varaut,{'INS'})~=0);
    genenow4=geneaut(isins);
    tum1now4=tum1aut(isins);
    tum2now4=tum2aut(isins);
    [unigene4,id,~] = unique(genenow4,'stable');
    for n4=1:length(id)
        isg3=find(ismember(genenow4,unigene4(n4))~=0);
        tum1now33=tum2now4{isg3};
        countc=length(find(tum1now33=='C'));
        countg=length(find(tum1now33=='G'));
        countt=length(find(tum1now33=='T'));
        counta=length(find(tum1now33=='A'));
        xi(n4)=(counta+countg)-(countc+countt);
        yi(n4)=(counta+countc)-(countg+countt);
        zi(n4)=(counta+countt)-(countg+countc); 
    end
    if n4~=0
        insscore=(xi.^2+yi.^2+zi.^2+1)/4;
        [~,loc4]=ismember(unigene4,geneall);
    else
        insscore=zeros(1,length(geneall));
        loc4=1:length(geneall);
    end
    geneins(i,loc4)=insscore;
    clear unigene4 insscore loc4 xi yi zi
end

save Gtotalscore.mat geneallmu genealltt genedel geneins geneall pid

```

> TTZ-feature=TMB+TTR+Z

```
clear
load Gtotalscore.mat geneallmu genealltt genedel geneins geneall pid

vscore=genedel+geneins;
[gmu,mu1,mu2]=rescaling(geneallmu,'minmax');
[gtt,tt1,tt2]=rescaling(genealltt,'minmax');
[gv,v1,v2]=rescaling(vscore,'minmax');
%save rescalingpara.mat mu1 mu2 tt1 tt2 v1 v2 geneall
delscore=genedel;
insscore=geneins;
scoreall=gmu+gtt+gv;
save gscoreanalysis.mat scoreall delscore insscore gmu gtt gv geneall pid vscore
iszero=find(all(scoreall==0,1)==0);%非全0列的Index
scorecandi=scoreall(:,iszero);
genecandi=geneall(iszero);
hscore=scorecandi((1:154),:);%16556个非全0列基因
nscore=scorecandi((155:end),:);
save gscore.mat hscore nscore scorecandi genecandi pid
```

- ```
  ######function-rescaling
  function [X,para1,para2]=rescaling(X,method,para1,para2)
  %+++   data pretreatment
  %+++ HD Li, Central South University
  
  
  if nargin==2
    [Mx,Nx]=size(X);
     if strcmp(method,'autoscaling')
      para1=mean(X);para2=std(X);
     elseif strcmp(method,'center')
      para1=mean(X);para2=ones(1,Nx);
     elseif strcmp(method,'unilength')
      para1=mean(X);
      for j=1:size(X,2);
      para2(1,j)=norm(X(:,j)-para1(j));
      end
     elseif strcmp(method,'minmax')
      para1=min(X);maxv=max(X);
      para2=maxv-para1;  
     elseif strcmp(method,'pareto');
      para1=mean(X);para2=sqrt(std(X));
     else
      display('Wrong data pretreat method!');
     end
     
     for i=1:Nx
        if para2(i)==0
            X(:,i)=zeros(Mx,1);
        else
       X(:,i)=(X(:,i)-para1(i))/para2(i);
        end
     end
     
  elseif nargin==4
     [Mx,Nx]=size(X);
     for i=1:Nx    
         if para2(i)==0
            X(:,i)=zeros(Mx,1);
        else
       X(:,i)=(X(:,i)-para1(i))/para2(i);
         end
     end
  end
  ```

- ```
  ######判断是否服从正态分布
  Hh=[];
  Hn=[];
  for v=1:16556
  A1 = hscore(:,v);
  A2 = nscore(:,v);
  alpha = 0.05;
  % 正态分布判断
  [mu1, sigma1] = normfit(A1);
  [mu2, sigma2] = normfit(A2);
  p1 = normcdf(A1, mu1, sigma1);
  p2 = normcdf(A2, mu2, sigma2);
  [H1,s1] = kstest(A1, [A1, p1], alpha);%H1=0 该数据源服从正态分布
  [H2,s2] = kstest(A2, [A2, p2], alpha);
  Hh=[Hh H1];
  Hn=[Hn H2];
  end
  ```

- ```
  %% Mann-Whitney U test ：
  %this test can be used to determine whether two independent samples were selected from populations 
  %having the same distribution;Unlike the t-test it does not require the assumption of normal distributions. 
  %It is nearly as efficient as the t-test on normal distributions
  pm=[];hm=[];
  for u=1:16556
      [pm2,hm2]=ranksum(hscore(:,u),nscore(:,u));
      pm=[pm pm2];
      hm=[hm hm2];
  end
  isone=find(hm==1);
  Xgcal=score(:,isone);
  ygcal=[ones(155,1);-ones(63,1)];
  sigggene=gene2(isone);
  ggeneall=geneall;
  save Gdatatotrain.mat gscoreall ggeneall sigggene
  ```

------

> 