function [MWUboots] = runmwuallboots(bootcontrols,experiment,tauniquenames,tauniqueindices)
%Run the bootstrapreplicates
[m,n]=size(bootcontrols);
MWUboots=zeros(length(tauniqueindices),n);
for i=1:n
    [MWUpvalvector] = MWUbytransposon(experiment,bootcontrols(:,i),tauniquenames,tauniqueindices);
    MWUboots(:,i)=MWUpvalvector;

end
end
