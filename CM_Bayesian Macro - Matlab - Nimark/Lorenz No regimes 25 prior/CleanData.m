function [Z,CPImean,NGDPmean,NonZeroCPI,NonZeroNGDP]=CleanData(maxsurv)

%clean up the data
% clear all
% clc
% maxsurv=25;


load('NGDPsurv');
load('NGDPactual');
load('CPIsurv');
load('dCPI');
load('logGDP');
load('dTFP');
load('FFunds');
dTFP=0.0025*dTFP;

% hpCPI=hpfilter(dCPI,1000);
trendCPI=dCPI-detrend(dCPI);
% ZFFunds=detrend(FFunds)*0.0025;
ZFFunds=FFunds-trendCPI(1:117);
ZFFunds=0.01*(ZFFunds-mean(ZFFunds));
hpGDP=hpfilter(logGDP,1400);
hpGDP=[hpGDP;hpGDP(end)+hpGDP(end)-hpGDP(end-1);];



for t=1:length(dTFP)-1;
    TFP(t)=sum(dTFP(1:t+1));
end
TFPdev=detrend(TFP);

TT=length(NGDPsurv);
SampleDates=[1981.5:.25:2010.75];

SampleNGDP=[SampleDates;zeros(50,length(SampleDates));];



SampleCPI=[SampleDates;zeros(50,length(SampleDates));];

dateold=1;
for j=1:TT;
    
    date=CPIsurv(j,1)+CPIsurv(j,2)*0.25 - 0.25;
    t=(date-1981.5)*4+1;
    
    if date == dateold;
        c=c+1 ;     
    else
        c=0;
    end
    
    
    if CPIsurv(j,7) ==-999;
        c=c-1;
    else
        SampleCPI(c+2,t)=0.01*CPIsurv(j,7)-0.01*trendCPI(t);
    end
     dateold=date;
end

CPImean=mean(mean(SampleCPI(2:14,:),2));

dateold=1;
for j=1:TT;
    
    date=NGDPsurv(j,1)+NGDPsurv(j,2)*0.25 - 0.25;
    t=(date-1981.5)*4+1;
    
    if date == dateold;
        c=c+1 ;     
    else
        c=0;
    end
    
    
    if NGDPsurv(j,7) ==-999 || NGDPsurv(j,6) ==-999 ;
        c=c-1;
    else
        SampleNGDP(c+2,t)=log(NGDPsurv(j,7)) - log(NGDPsurv(j,6 ))-0.0025*trendCPI(t)-0.0025*CPImean-(hpGDP(t+1)-hpGDP(t));
    end
     dateold=date;
end

NGDPmean=0.0234;%mean(mean(SampleNGDP(2:7,:),2));


NonZeroNGDP=zeros(length(SampleNGDP),1);
for t=1:length(SampleNGDP);
    NonZeroNGDP(t)=sum(SampleNGDP(2:end,t)~=0);
    NonZeroNGDP(t)=min(maxsurv,NonZeroNGDP(t));
end

NonZeroCPI=zeros(length(SampleCPI),1);
for t=1:length(SampleNGDP);
    NonZeroCPI(t)=sum(SampleCPI(2:end,t)~=0);
    NonZeroCPI(t)=min(maxsurv,NonZeroCPI(t));
end
% for t=1:length(SampleNGDP);
%     hist(SampleNGDP(2:1+NonZeroNGDP(t),t),25);
%     t
%     pause
% end
% 
% for t=1:length(SampleCPI);
%     hist(SampleCPI(2:1+NonZeroCPI(t),t),25);
%     t
%     pause
% end
% 


ZTFP=TFP'-hpfilter(TFP,1400);
ZCPI=detrend(0.01*dCPI(1:117));
ZGDP=logGDP(1:117)-hpGDP(1:117);
ZCPIsurv=SampleCPI(2:end,1:117);
ZGDPsurv=SampleNGDP(2:end,1:117);
Z=[ZTFP';0.25*ZFFunds';0.25*ZCPI';ZGDP';0.25*(ZCPIsurv(1:maxsurv,:));(ZGDPsurv(1:maxsurv,:));];
% 
% figure(1)
% plot(Z(1:4,:)')
% 
% figure(2)
% plot(Z(3,:)','linewidth',3)
% hold on
% plot(0.25*(ZCPIsurv(1:maxsurv,:))')
% 
% figure(3)
% plot(Z(3,:)','linewidth',3)
% hold on
% plot(Z(4,:)','linewidth',3,'linestyle',':')
% hold on
% plot((ZGDPsurv(1:maxsurv,:))')
% 
% S=cov(Z(1:4,:)')
% r=corrvc(S)

    
