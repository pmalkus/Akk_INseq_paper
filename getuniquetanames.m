function [uniquenames,uniqueindices]=getuniquetanames(names)
j=1;
h=1;
for i=1:length(names)-1;
    match=strmatch(names(i,1),names(i+1));
   if isempty(match)==0;
       h=h+1;
   end
   if isempty(match)==1
       uniquenames(j,1)=names(i,1);
       uniqueindices(j,1)=(i-h)+1;
       uniqueindices(j,2)=i;
       j=j+1;
       h=1;
   end
   if i==length(names)-1
       uniquenames(j,1)=names(i,1);
       uniqueindices(j,1)=((i+1)-h)+1;
       uniqueindices(j,2)=i+1;
   end
end

       