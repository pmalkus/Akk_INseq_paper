function [MWUpvalvector] = MWUbytransposon(experiment,control,tauniquenames,tauniqueindices)
%
%   Detailed explanation goes here
for i=1:length(tauniquenames)
    labels=[];
    scores1=[];
    scores2=[];
    labels(1:((tauniqueindices(i,2)-tauniqueindices(i,1))+1),1)=0;
    labels(((tauniqueindices(i,2)-tauniqueindices(i,1))+2):(2*((tauniqueindices(i,2)-tauniqueindices(i,1))+1)),1)=1;
    scores1(1:((tauniqueindices(i,2)-tauniqueindices(i,1))+1),1)=experiment(tauniqueindices(i,1):tauniqueindices(i,2),1);
    scores2(1:((tauniqueindices(i,2)-tauniqueindices(i,1))+1),1)=control(tauniqueindices(i,1):tauniqueindices(i,2),1);
    [p,h] = ranksum(scores1,scores2);
    MWUpvalvector(i,1)=p;
end

