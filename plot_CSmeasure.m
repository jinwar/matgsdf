figure(1)
clf
hold on
for ip=1:length(periods)
for ics = 1:length(CS_measure)
  x(ics) = CS_measure(ics).ddist;
  y(ics) = CS_measure(ics).fitpara(4,ip);
end
para = polyfit(x,y,1);
title(['V: ',num2str(1./para(1)),' Offset: ',num2str(para(2))]);
plot(x,y,'x')
phaseV(ip)=1./para(1);
end
figure(2)
clf
plot(periods,phaseV,'x-')