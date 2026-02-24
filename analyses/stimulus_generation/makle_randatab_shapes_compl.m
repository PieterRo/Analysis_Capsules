close all
clear all
addpath(genpath('\\vs03\vs03-vandc-2\Object_attention\Shapes\_code\code_utils'));

monkey = 'monkeyF';
dir_gen = '\\vs03\VS03-VandC-2\Object_attention\Shapes\_stims\_old\';
dir_stims = ['\\vs03\VS03-VandC-2\Object_attention\Shapes\_stims\',monkey,'\'];
RFs_array = ['\\vs03\vs03-vandc-2\Object_attention\Shapes\_stims\',monkey,'\RFs\'];
logs_dir = ['\\vs03\vs03-vandc-2\Object_attention\Shapes\',monkey,'\_logs\'];

th = 6.5;
gc_conditions = [2,4,6];
classes = 2;
RF_to_cone = 0.75;
% arrays = [2,3,4,5,8,9,10,101];  %monkeyN
% arrays_notsel = 3; %monkeyN
arrays = [2,5,7]; %monkeyF
arrays_notsel = arrays; %monkeyF
HX = 1024/2;
HY = 768/2;
dist_2nd_cues = 14*25;
rgb{1} = [.6 .6 .7];
rgb{2} = [.6 .62 .47];
theta = 0:pi/180:2*pi;
angles = 0:.01:pi; 
complexities = 4;
nreps = 12;
make_curves = 0;
make_stims = 1;
makerand = 1;
iter_max = 500;
ratio_ = 38.05/44.82; %monkey F
% ratio_ = 1; %monkey N

if make_curves
    clear todo
    for class_ = 1:classes
        clear fold sp file image_ dir_ start_r start_c

        if class_ == 1
            fold = 'Animals\';
            sp = [0 2 10 100]; %crocodile
            [image_,dir_,start_r,start_c] = curv_img_proc([dir_gen,fold,'crocodile.png']);
        else
            fold = 'Animals\';
            sp = [0 3 15 100]; %monkey
            [image_,dir_,start_r,start_c] = curv_img_proc([dir_gen,fold,'monkey2.png']);
        end
        if ~isempty(dir_)
            clear chain cc
            chain = bwtraceboundary(image_,[start_r,start_c],dir_);
            [cc] = chaincode(chain);
            
            ind = 0;
            for j = sp
                ind = ind+1;
                if ind == 1
                    clear image_0 dir_0 start_r0 start_c0
                    %                     [image_0,dir_0,start_r0,start_c0] = curv_img_proc([dir_gen,fold,file{1}]);
                    [image_0,dir_0,start_r0,start_c0] = curv_img_proc([dir_gen,fold,'0.png']);
                    if ~isempty(dir_0)
                        clear chain0 cc0 x_ x
                        chain0 = bwtraceboundary(image_0,[start_r0,start_c0],dir_0);
                        [cc0] = chaincode(chain0);
                        x_ = fourier_approx(transpose([cc0.code]),100,4096,1);
                        x = [x_; x_(1,1) x_(1,2)];
                        todo(class_,ind,:,:) = x;
                    end
                else
                    x_ = fourier_approx(transpose([cc.code]),j,4096,1);
                    x = [x_; x_(1,1) x_(1,2)];
                    todo(class_,ind,:,:) = x;
                end
            end
            
            for p = 1:length(sp)
                clear name_z X Y X2d Y2d IN
                name_z = ['class',mat2str(class_),'_complexity',mat2str(p)];
                X = squeeze(todo(class_,p,:,1))*200+300;
                Y = squeeze(todo(class_,p,:,2))*200+150;
                [X2d,Y2d] = meshgrid(1:600,1:300);
                IN = inpolygon(X2d,Y2d,X,Y);
                
                figure
                clear all_centers all_radii X Y s idx
                [all_centers, all_radii, X, Y] = growth_cones_filling(IN,th,0);
                plot(X,Y,'k','LineWidth',2)
                hold on
                viscircles(all_centers,all_radii);
                set(gca,'DataAspectRatio', [1 1 1]);
                axis off
                title('Pick the 1st cue, then the 2nd cue and finally 3 cones:');
                drawnow
                
                [s(:,1),s(:,2)] = getpts;
                close all
                
                for ss = 1:length(s)
                    clear temp which_s
                    temp = s(ss,:);
                    which_s = find(sqrt((temp(1)-all_centers(:,1)).^2+(temp(2)-all_centers(:,2)).^2)<all_radii');
                    s(ss,:) = all_centers(which_s,:);
                    idx(ss) = which_s;
                end
                plot(X,Y,'k','LineWidth',2)
                hold on
                viscircles(all_centers,all_radii);
                scatter(s(1:2,1),s(1:2,2),'filled','b')
                scatter(s(3:end,1),s(3:end,2),'filled','g')
                set(gca,'DataAspectRatio', [1 1 1]);
                axis off
                
                clear cues RFs mean_size
                cues = s(1:2,:);
                RFs = s(3:end,:);
                mean_size = min(mean(all_radii(idx(3:end)))*RF_to_cone,min(all_radii(idx(3:end))));
                viscircles(all_centers(idx(3:end),:),repmat(mean_size,1,length(idx(3:end))),'Color','g');
                drawnow
                temp_figures.(name_z).p = p;
                temp_figures.(name_z).class = class_;
                temp_figures.(name_z).X = X;
                temp_figures.(name_z).Y = Y;
                temp_figures.(name_z).cues = cues;
                temp_figures.(name_z).RFs = RFs;
                temp_figures.(name_z).idx = idx;
                temp_figures.(name_z).all_centers = all_centers;
                temp_figures.(name_z).all_radii = all_radii;
                temp_figures.(name_z).mean_size = mean_size;
                
            end
            set(gca,'DataAspectRatio', [1 1 1]);
        end
    end
else
    load([logs_dir,'temp_figs.mat']);
end

if make_stims
    for p = 1:complexities
%     for p = 4
        disp(['%%% %%% Complexity ',num2str(p)])
        clear ALLMAT ALLCOORDS RANDTAB
        stim_dir = [dir_stims,'comp',mat2str(p),'\'];
        if ~exist(stim_dir, 'dir')
            mkdir(stim_dir);
        end
        exp_stim_dir = [dir_stims,'exp\comp',mat2str(p),'\'];
        if ~exist(exp_stim_dir, 'dir')
            mkdir(exp_stim_dir);
        end
%         if p < 4
%             angles = 0:.01:pi/2; %comp 2-3
%         else
%             angles = 0:.01:pi; %comp 4
%         end
        z = 1;
        for array = arrays
            disp(['%%% Array ',num2str(array)])
            clear temp_cents array_RF_x array_RF_y array_RF_sz best_sizes
            if ~ismember(array,arrays_notsel)
                load([RFs_array,'array',num2str(array),'_good_rfs_selection.mat'])
                temp_cents = ratio_.*nanmedian([centrex(fin_goods);centrey(fin_goods)],2);
                array_RF_x  = ratio_.*temp_cents(1);
                array_RF_y  = ratio_.*temp_cents(2);
                array_RF_radius  = ratio_.*(median(mean([szx(fin_goods);szy(fin_goods)]))/2);
            else
                load([RFs_array,'array',num2str(array),'_good_rfs.mat'])
                temp_cents = ratio_.*nanmedian([centrex;centrey],2);
                array_RF_x  = ratio_.*temp_cents(1);
                array_RF_y  = ratio_.*temp_cents(2);
                array_RF_radius  = ratio_.*(median(mean([szx;szy]))/2);
            end
            name_array = ['array_',num2str(array)];
            for class_ = 1:classes
                clear name_z X Y cues RFs idx all_centers all_radii mean_size
                name_z = ['class',mat2str(class_),'_complexity',mat2str(p)];
                X = temp_figures.(name_z).X;
                Y = temp_figures.(name_z).Y;
                cues = temp_figures.(name_z).cues;
                RFs = temp_figures.(name_z).RFs;
                idx = temp_figures.(name_z).idx;
                all_centers = temp_figures.(name_z).all_centers;
                all_radii = temp_figures.(name_z).all_radii;
                mean_size = temp_figures.(name_z).mean_size;
                
                for cc = 1:length(cues)
                    clear s t
                    s = cues(cc,:);
                    t = cues(3-cc,:);
                    
                    for rr = 1:length(RFs)
                        cost_bads2 = 1;
                        iter = 0;
                        while cost_bads2==1
                            clear cones temp_ temp_ratio temp_center temp_X temp_Y
                            if cc == 1
                                cones = gc_conditions(rr);
                            else
                                cones = gc_conditions(4-rr);
                            end
                            name_c = ['cones',mat2str(cones)];
                            temp_ = RFs(rr,:);
                            
                            temp_X = [X; s(1); t(1); temp_(1)];
                            temp_Y = [Y; s(2); t(2); temp_(2)];
                            angle = angles(randperm(length(angles),1));
                            
                            R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
                            A = R*[temp_X,temp_Y]';
                            temp_X = A(1,:)';
                            temp_Y = A(2,:)';
                            temp_ = [temp_X(end) temp_Y(end)];
                            temp_X(end) = [];
                            temp_Y(end) = [];
                            temp_t = [temp_X(end) temp_Y(end)];
                            temp_X(end) = [];
                            temp_Y(end) = [];
                            temp_s = [temp_X(end) temp_Y(end)];
                            temp_X(end) = [];
                            temp_Y(end) = [];
                            
                            temp_ratio = array_RF_radius/mean_size;
                            temp_center = temp_.*temp_ratio;
                            temp_X = temp_X.*temp_ratio - (temp_center(1) - array_RF_x);
                            temp_Y = temp_Y.*temp_ratio - temp_center(2) + array_RF_y;
                            cost_bads = sum(~inpoly([temp_X temp_Y],[-HX -HY; HX -HY; HX HY; -HX HY]));
                            iter = iter+1;
                            if ~cost_bads
                                temp_s = temp_s.*temp_ratio - (temp_center - [array_RF_x array_RF_y]);
                                temp_t = temp_t.*temp_ratio - (temp_center - [array_RF_x array_RF_y]);
                                
                                dist_cue_temp = norm(temp_s-temp_t);
                                [a,b] = circcirc(temp_s(1),temp_s(2),dist_cue_temp,temp_t(1),temp_t(2),dist_2nd_cues);
                                idx_2nd = randperm(2,1);
                                temp_t2 = [a(idx_2nd) b(idx_2nd)];
                                
                                clear name_z2 X2 Y2 cues2
                                name_z2 = ['class',mat2str(3-class_),'_complexity',mat2str(p)];
                                X2 = temp_figures.(name_z2).X;
                                Y2 = temp_figures.(name_z2).Y;
                                cues2 = temp_figures.(name_z2).cues;
                                
                                xcircle_RF = array_RF_radius * cos(theta) + array_RF_x;
                                ycircle_RF = array_RF_radius * sin(theta) + array_RF_y;
                                
                                for ref = 1:2
                                    [~,~,transform] = procrustes([temp_s;temp_t2],cues2,'reflection',ref-1);
                                    Z_points = transform.b*[X2 Y2]*transform.T + transform.c(1,:);
                                    temp_X2 = Z_points(:,1);
                                    temp_Y2 = Z_points(:,2);
                                    cost_RF(ref) = sum(inpoly([xcircle_RF;ycircle_RF]',[temp_X2,temp_Y2]));
                                end
                                ref = find(cost_RF==min(cost_RF),1);
                                [~,~,transform] = procrustes([temp_s;temp_t2],cues2,'reflection',ref-1);
                                Z_points = transform.b*[X2 Y2]*transform.T + transform.c(1,:);
                                temp_X2 = Z_points(:,1);
                                temp_Y2 = Z_points(:,2);
                                cost_bads2 = min(cost_RF) + sum(~inpoly([temp_X2 temp_Y2],[-HX -HY; HX -HY; HX HY; -HX HY]));
                                
                                
                                if ~cost_bads2
                                    angle = rad2deg(angle);
                                    for fg = 1:2
                                        for col = 1:2
                                            clear col_fig col_back
                                            if col == 1
                                                col_fig = rgb{1};
                                                col_back = rgb{2};
                                            else
                                                col_fig = rgb{2};
                                                col_back = rgb{1};
                                            end
                                            
                                            clear temp_t_dis temp_t_fig X_fig Y_fig X_back Y_back
                                            if fg == 1 %fg == 1 (back on top) fg == 2 (array on top)
                                                temp_t_dis = temp_t;
                                                temp_t_fig = temp_t2;
                                                X_fig = temp_X2;
                                                Y_fig = temp_Y2;
                                                X_back = temp_X;
                                                Y_back = temp_Y;
                                            else
                                                temp_t_dis = temp_t2;
                                                temp_t_fig = temp_t;
                                                X_fig = temp_X;
                                                Y_fig = temp_Y;
                                                X_back = temp_X2;
                                                Y_back = temp_Y2;
                                            end
                                            
                                            clear bw_fig bw_back fig_temp back_temp final
                                            bw_fig = poly2mask(double(X_fig)+HX,HY-double(Y_fig),HY*2,HX*2);
                                            bw_back = poly2mask(double(X_back)+HX,HY-double(Y_back),HY*2,HX*2);
                                            bw_back(bw_fig) = 0;
                                            fig_temp = double(cat(3, bw_fig.*col_fig(1), bw_fig.*col_fig(2), bw_fig.*col_fig(3)));
                                            back_temp = double(cat(3, bw_back.*col_back(1), bw_back.*col_back(2), bw_back.*col_back(3)));
                                            final = fig_temp + back_temp;
                                            final(final == 0) = .5;
                                            
                                            imshow(final)
                                            hold on
                                            scatter(HX,HY,'filled','r')
                                            plot(xcircle_RF+HX,HY-ycircle_RF,'g','LineWidth',2);
                                            scatter(temp_s(1)+HX,HY-temp_s(2),'filled','y')
                                            scatter(temp_t(1)+HX,HY-temp_t(2),'filled','y')
                                            scatter(temp_t2(1)+HX,HY-temp_t2(2),'filled','y')
                                            set(gca,'DataAspectRatio', [1 1 1]);
                                            drawnow
                                            
                                            
                                            clear name_z2
                                            name_z2 = ['stim_',num2str(z)];
                                            saveas(gcf,[exp_stim_dir,name_array,'_',name_z,'_',name_c,'_',name_z2,'.png'])
                                            pause(.5)
                                            close all
                                            
                                            ALLMAT(z,:) = [z array cones fg class_ angle cc col];
                                            ALLCOORDS.(name_z2).s = temp_s;
                                            ALLCOORDS.(name_z2).t_fig = temp_t_fig;
                                            ALLCOORDS.(name_z2).t_back = temp_t_dis;
                                            imwrite(final,[stim_dir,'\',sprintf('%03.0f.bmp', z)]);
                                            pause(.5)
                                            z = z+1;
                                        end
                                    end
                                else
                                    cost_bads2 =1;
                                    if iter >= iter_max
                                        cost_bads2 = 2;
                                        disp('Max iter 2!')
                                    end
                                    %                                     disp('bad #2!')
                                end
                            else
                                cost_bads2 =1;
                                if iter >= iter_max
                                    cost_bads2 = 2;
                                    disp('Max iter 1!')
                                    cones
                                end
%                                 disp('bad #1!')
                            end
                        end
                    end
                end
            end
        end
        
        ALLMAT_reps = [];
        for i = 1:nreps
            ALLMAT_reps = [ALLMAT_reps; ALLMAT];
        end
        if makerand
            perm_idx = randperm(length(ALLMAT_reps),length(ALLMAT_reps));
            RANDTAB = ALLMAT_reps(perm_idx,:);
        else
            RANDTAB = ALLMAT_reps;
        end
        
        filename = [logs_dir,'RANDTAB_shapes_comp',mat2str(p),'_',monkey,'.mat'];
        filename_back = [logs_dir,'RANDTAB_shapes_comp',mat2str(p),'_',monkey,'_backup.mat'];
        
        save(filename,'ALLMAT','RANDTAB','ALLCOORDS','ALLMAT_reps','perm_idx');
        pause(.5)
        save(filename_back,'ALLMAT','RANDTAB','ALLCOORDS','ALLMAT_reps','perm_idx');
        pause(.5)
        
        sum(ALLMAT(:,3)==2)
        sum(ALLMAT(:,3)==6)
        sum(ALLMAT(:,3)==4)
    end
end

