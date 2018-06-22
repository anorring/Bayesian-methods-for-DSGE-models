%Program to compute correlations between changes and dispersion
clear all
clc

load CPId
load CPIdisp
load NGDPd
load NGDPdisp

samplecut=50;

CPId=CPId(samplecut:end);
CPIdisp=CPIdisp(samplecut:end);
NGDPd=NGDPd(samplecut:end);
NGDPdisp=NGDPdisp(samplecut:end);

%normalize
CPId=CPId-mean(CPId);
CPIdisp=CPIdisp-mean(CPIdisp);
NGDPd=NGDPd-mean(NGDPd);
NGDPdisp=NGDPdisp-mean(NGDPdisp);

CPIdabs=abs(CPId)-mean(abs(CPId));
NGDPdabs=abs(NGDPd)-mean(abs(NGDPd));

CPId=CPId./(cov(CPId)^0.5);
CPIdisp=CPIdisp./(cov(CPIdisp)^0.5);


NGDPd=NGDPd./(cov(NGDPd)^0.5);
NGDPdisp=NGDPdisp./(cov(NGDPdisp)^0.5);

CPIdabs=CPIdabs./(cov(CPIdabs)^0.5);

NGDPdabs=NGDPdabs./(cov(NGDPdabs)^0.5);

infld=CPId(2:end)-CPId(1:end-1);
infld=infld-mean(infld);

inflabs=abs(infld)-mean(abs(infld));
inflabs=inflabs./(cov(inflabs)^0.5);
infld=infld./(cov(infld)^0.5);

subplot(2,1,1)
plot(CPId);hold on; plot(CPIdisp,'k')

subplot(2,1,2)
plot(NGDPd);hold on; plot(NGDPdisp,'k')

cov(CPId,CPIdisp)
cov(CPIdabs,CPIdisp)

cov(NGDPd,NGDPdisp)
cov(NGDPdabs,NGDPdisp)

cov(infld,CPIdisp(2:end))
cov(inflabs,CPIdisp(2:end))



