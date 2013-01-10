function winpara = auto_win_select(event,periods,mingroupv,maxgroupv,bandnum,center_freq)
% This function is used to automatically select the window range used for gsdf method.
% The output format is
% v1 = winpara(1); t1 = winpara(2);
% v2 = winpara(3); t2 = winpara(4);
% and the window is defined by L/v1+t1 -- L/v2+t2

if ~exist('mingroupv')
    mingroupv = 2;
end
if ~exist('maxgroupv')
    maxgroupv = 5;
end
if ~exist('bandnum')
    bandnum = 20;
end
if ~exist('center_freq')
    center_freq = 0.025;
end

cycle_before = 2;
cycle_after = 5;

peakamptol = 0.5;   % the ratio of accepted peak compare to the largest peak
peak_search_range = 1;
groupv_diff_tol = 0.15;
positive_disp_weight = 5;


minf = 1/periods(end);
maxf = 1/periods(1);
freqs = linspace(minf,maxf,bandnum);
[temp center_freq_index] = min(abs(freqs - center_freq));

isdebug =1;

for ista = 1:length(event.stadata)
%  for ista = 20
    ista
    % set up time axis
    bgtime = event.stadata(ista).otime - event.otime;
    dt = event.stadata(ista).delta;
    Nt = length(event.stadata(ista).data);
    taxis = bgtime + [0:Nt-1]'*dt;
    
    % Build up gaussian filter
    [gausf,faxis] = build_gaus_filter(freqs,dt,Nt,0.06,0.1);
    
    % get original data and make the fourier transform
    odata = event.stadata(ista).data;
    if size(odata,1) == 1  % in matlab, the fast direction is column
        odata = odata';
    end
    fftodata = fft(odata);
    
    clear envelop_nbands nband nbands norm_envelop
    % apply narrow-band filters
    for ip = 1:length(freqs);
        nband = fftodata .* [gausf(:,ip); zeros(Nt-length(gausf(:,ip)),1)];
        nband = ifft(nband);
        nbands(:,ip) = nband;
    end % end of loop ip
    
    % remove the area out of window defined by min and max group velocity
    tmin = event.stadata(ista).dist/maxgroupv;
    tmax = event.stadata(ista).dist/mingroupv;
    if tmin < taxis(1) || tmax > taxis(end)
        disp(['Station ',event.stadata(ista).stnm,' does not contain enough data']);
        continue;
    end
    inwin_ind = find(taxis > tmin & taxis < tmax);
    outwin_ind = find(taxis < tmin | taxis > tmax);
    envelop_nbands=abs(nbands);
    envelop_nbands(outwin_ind,:) = 0; % apply group velocity window
    for ip = 1:length(freqs) % Normalize it for grading purpose
        envelop_nbands(:,ip) = envelop_nbands(:,ip) / max(envelop_nbands(:,ip));
    end
    
    
    % Find local maximum as potential peaks
    clear peaks
    for ip = 1:length(freqs)
        envelopfun = envelop_nbands(:,ip);
        diff_fun = diff(envelopfun);
        peakind = inwin_ind(find(diff_fun(inwin_ind-1)>0 & diff_fun(inwin_ind)<0));
        bigpeakamp = max(envelopfun);
        peakind = peakind(find(envelopfun(peakind)>bigpeakamp*peakamptol));
        peaks(ip).peaktimes = taxis(peakind);
        peaks(ip).peakamps = envelopfun(peakind);
        peaks(ip).peaknum = length(peakind);
    end
    
    % From center frequency band, search for the best peak to both sides
    ip = center_freq_index;
    [temp ind] = max(peaks(ip).peakamps);
    peaks(ip).bestpeak = peaks(ip).peaktimes(ind);
    peaks(ip).bestpeaksnr = peaks(ip).peakamps(ind)/sum([peaks(ip).peakamps]);
    % To the lower frequency bands
    for ip=center_freq_index-1:-1:1  % loop for lower frequencies
        [temp closest_peak_i] = min(abs(peaks(ip+1).bestpeak - peaks(ip).peaktimes));
        bestpeaki = closest_peak_i;
        bestpeakamp = peaks(ip).peakamps(bestpeaki);
        groupV0 = event.stadata(ista).dist / peaks(ip+1).bestpeak;
        for ipeak = closest_peak_i-peak_search_range:closest_peak_i+peak_search_range
            if ipeak >= 1 && ipeak <= length(peaks(ip).peaktimes) % search valid
                groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(ipeak);
                if groupvi > groupV0*(1-groupv_diff_tol) && ...
                        groupvi < groupV0*(1+groupv_diff_tol)
                    if peaks(ip).peakamps(ipeak) > bestpeakamp
                        bestpeaki = ipeak;
                        bestpeakamp = peaks(ip).peakamps(ipeak);
                    end % end for finding largest nearby peak
                end % end for check group v
            end % end for check ipeak vavid
        end % end of loop ipeak
        peaks(ip).bestpeak = peaks(ip).peaktimes(bestpeaki);
        peaks(ip).bestpeaksnr = peaks(ip).peakamps(bestpeaki)/sum([peaks(ip).peakamps]);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        if groupvi < groupV0*(1-groupv_diff_tol) || ...
                groupvi > groupV0*(1+groupv_diff_tol)
            peaks(ip).bestpeak = peaks(ip+1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
        end
    end % end of loop ip
    % To the higher frequency bands
    for ip=center_freq_index+1:length(freqs)  % loop for higher frequencies
        [temp closest_peak_i] = min(abs(peaks(ip-1).bestpeak - peaks(ip).peaktimes));
        bestpeaki = closest_peak_i;
        amp_point = peaks(ip).peakamps(bestpeaki);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        groupV0 = event.stadata(ista).dist / peaks(ip-1).bestpeak;
        if groupvi/groupV0 > 1
            v_point = 1-(groupvi/groupV0-1)/groupv_diff_tol;
        else
            % double the weight for having a positive dispersion curve
            v_point = 1-(1-groupvi/groupV0)/groupv_diff_tol/positive_disp_weight;
        end
        best_sum_point = amp_point + v_point;
        for ipeak = closest_peak_i-peak_search_range:closest_peak_i+peak_search_range
            if ipeak >= 1 && ipeak <= length(peaks(ip).peaktimes) % search valid
                groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(ipeak);
                if groupvi > groupV0*(1-groupv_diff_tol) && ...
                        groupvi < groupV0*(1+groupv_diff_tol)
                    amp_point = peaks(ip).peakamps(ipeak);
                    if groupvi/groupV0 > 1
                        v_point = 1-(groupvi/groupV0-1)*10;
                    else
                        % double the weight for having a positive dispersion curve
                        v_point = 1-(1-groupvi/groupV0)*10/positive_disp_weight;
                    end
%                     disp([num2str(ip),':',num2str(amp_point),',',num2str(v_point)]);
                    if amp_point + v_point > best_sum_point
                        bestpeaki = ipeak;
                        best_sum_point = amp_point + v_point;
                    end % end for finding largest nearby peak
                end % end for check group v
            end % end for check ipeak vavid
        end % end of loop ipeak
        peaks(ip).bestpeak = peaks(ip).peaktimes(bestpeaki);
        peaks(ip).bestpeaksnr = peaks(ip).peakamps(bestpeaki)/sum([peaks(ip).peakamps]);
        groupvi = event.stadata(ista).dist / peaks(ip).peaktimes(bestpeaki);
        if groupvi < groupV0*(1-groupv_diff_tol) || ...
                groupvi > groupV0*(1+groupv_diff_tol)
            peaks(ip).bestpeak = peaks(ip-1).bestpeak;
            peaks(ip).bestpeaksnr = 0;
        end
    end % end of loop ip
    
    if 0
        for ip = 1:length(freqs)
            norm_envelop(:,ip) = envelop_nbands(:,ip) / max(envelop_nbands(:,ip));
        end
        figure(36)
        clf
        hold on
        [xi yi] = ndgrid(taxis,freqs);
        contourf(xi,yi,norm_envelop);
        for ip = 1:length(freqs)
            plot(peaks(ip).bestpeak,freqs(ip),'gx','markersize',15,'linewidth',3);
        end
        for ip = 1:length(freqs)
            plot(peaks(ip).peaktimes,ones(size(peaks(ip).peaktimes))*freqs(ip),'r.','markersize',15);
        end
    end
    
    % record the group delay
    groupdelay(:,ista) = [peaks(:).bestpeak]';
    snr(:,ista) = [peaks(:).bestpeaksnr]';
end % end of loop sta
dist = [event.stadata(:).dist];
for ip = 1:length(freqs)
    ip
    para = wpolyfit(dist,groupdelay(ip,:),snr(ip,:),1);
    groupv(ip) = 1/para(1);
    offset(ip) = para(2);
    if isdebug
        figure(38)
        clf
        hold on
        plot(dist,groupdelay(ip,:),'x');
        plot([min(dist),max(dist)],...
            [min(dist)/groupv(ip)+offset(ip),max(dist)/groupv(ip)+offset(ip)])
        title(['V:',num2str(groupv(ip)),' t:',num2str(offset(ip))]);
        pause
    end
end % end of ip
winpara = groupv;
end % end of function
