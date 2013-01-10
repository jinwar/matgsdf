function [groupv offset] = groupv_fit(dist,time,snr,mingroupv,maxgroupv)
% function to apply weighted least square polynomial fitting to the data
%
x=dist(:);
y=time(:);
w = snr(:);
N=1;

if ~exist('mingroupv')
    mingroupv = 2;
end
if ~exist('maxgroupv')
    maxgroupv = 5;
end

if length(x)~=length(y) || length(x)~=length(w)
	disp('The input vectors should have same length!');
	return 
end
mat = zeros(length(x),N+1);
for i=1:length(x)
	for j=1:N+1
		mat(i,j) = x(i)^(N+1-j);
	end
end
W = diag(w);

para = (mat'*W*mat)\(mat'*W*y);

groupv = 1/para(1);
offset = para(2);

if groupv > maxgroupv
	groupv = maxgroupv;
	offset = sum((y - x/groupv).*w)/sum(w);
elseif groupv < mingroupv;
	groupv = mingroupv;
	offset = sum((y - x/groupv).*w)/sum(w);
end


end


