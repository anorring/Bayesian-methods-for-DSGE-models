%transform Doms and Morin data
clear all

load('DMgoodraw');
load('DMbadraw');

sampleDates=[1981.5:.25:2003.6];


DMgood=[];
for j=1981.5:.25:2010.5;
    s=[];
    
    for t=1:length(DMgoodraw);
    if DMgoodraw(t,1) >= j && DMgoodraw(t,1) < j+0.25;
        s=[s DMgoodraw(t,2)];
    
    end
    end
    DMgood=[DMgood mean(s)];
end
        

DMbad=[];
for j=1981.5:.25:2010.5;
    s=[];
    
    for t=1:length(DMbadraw);
    if DMbadraw(t,1) >= j && DMbadraw(t,1) < j+0.25;
        s=[s DMbadraw(t,2)];
    
    end
    end
    DMbad=[DMbad mean(s)];
end


DMgood=DMgood-min(DMgood);
DMgood=DMgood/max(DMgood);

DMbad=DMbad-min(DMbad);
DMbad=DMbad/max(DMbad);

plot(DMgood);
hold on
plot(DMbad,'k');

save('DMgood','DMgood');
save('DMbad','DMbad');
