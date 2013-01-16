
for ip=1:length(periods)
figure(1)
clf
hold on
for ics = 1:length(CS)
  x(ics) = CS(ics).ddist;
  y(ics) = CS(ics).dtp(ip);
  isgood(ics) = CS(ics).isgood(ip);
end
goodind = find(isgood>0);
para = polyfit(x(goodind),y(goodind),1);
title(['Period:',num2str(periods(ip)),' V: ',num2str(1./para(1)),' Offset: ',num2str(para(2))]);
plot(x(goodind),y(goodind),'x')
badind = find(isgood==ErrorCode.low_cohere);
plot(x(badind),y(badind),'rx')
badind = find(isgood==ErrorCode.high_tp_err);
plot(x(badind),y(badind),'gx')
phaseV(ip)=1./para(1);
plot([min(x) max(x)],[min(x)*para(1)+para(2), max(x)*para(1)+para(2)],'r');
pause
end
figure(2)
clf
plot(periods,phaseV,'x-')