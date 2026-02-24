% 
% %v1 signal
% for gc = 1:length(gc_cones)
%     for fg = 1:2
%         temp_signal = [];
%         for roi = 1:length(rois_)
%             if rois_(roi) <= 8
%                 clear f mean_temp
%                 roi_chns = array_chns{rois_(roi)};
%                 roi_chns = (floor(roi_chns(1)/64)*64+1):(1+floor(roi_chns(1)/64))*64;
%                 snr_chns = SNR_day(roi_chns) >= snr_th;
%                 roi_chns = roi_chns(snr_chns);
%                 f = ALLMAT(:,3) == gc_cones(gc) & corrs & ALLMAT(:,2) == rois_(roi) & ALLMAT(:,4) == fg;
%                 mean_temp = squeeze(nanmean(nanmean(squeeze(normMUA(roi_chns,f,:)))));
%                 if fg == 1 && ~isnan(sum(mean_temp))
%                     for chn = roi_chns
%                         clear temp_ground temp_figure temp_att temp_att_pad lat rs
%                         f = ALLMAT(:,3) == gc_cones(gc) & corrs & ALLMAT(:,2) == rois_(roi) & ALLMAT(:,4) == 1;
%                         temp_ground = squeeze(nanmean(squeeze(normMUA(chn,f,:))));
%                         f = ALLMAT(:,3) == gc_cones(gc) & corrs & ALLMAT(:,2) == rois_(roi) & ALLMAT(:,4) == 2;
%                         temp_figure = squeeze(nanmean(squeeze(normMUA(chn,f,:))));
%                         temp_att = temp_figure-temp_ground;
%                             temp_att_pad = padarray(temp_att,[0 200],nanmean(temp_att(gt_long)),'post');
%                             [lat,~,rs] = latencyfit4AM(smooth(temp_att_pad-mean(temp_att_pad(1:200)),50),tb_pad,2,0);
%                             if rs>.5
%                                 all_lats(gc,chn) = lat;
%                             end
%                             all_resps(gc,chn) = nanmean(temp_att(gt));
%                     end
%                 end
%                 temp_signal = [temp_signal mean_temp];
%             end
%         end
%         v1_signal(gc,fg,:) = smooth(nanmean(temp_signal,2),20);
%     end
% end
