function step1_setup_analysis(data_path, sample_num, skip, is_test)
    temp = dir(fullfile(data_path,'v*'));
    View_name = {temp.name};
    view_num = size(View_name, 2);
    for i = 1:view_num
        view_path = fullfile(data_path,View_name{i});
        if ~is_test && skip==0
            delete(fullfile(view_path,'*.csv'));
        end
        temp = dir(fullfile(view_path,'MIP_ROI_*'));
        Mip_name =  natsortfiles({temp.name});
        mip_num = size(Mip_name, 2);
        for j = skip+1:mip_num
            mip_path = fullfile(view_path,Mip_name{j});
            temp = dir(fullfile(mip_path,'ROI*'));
            Roi_name = natsortfiles({temp.name});
            roi_num = size(Roi_name, 2);
            for k = 1:roi_num
                roi_path = fullfile(mip_path,Roi_name{k});
                fprintf(datestr(now()));
                fprintf([' ', strrep(roi_path,'\','/'),'\n']);
                
                flag = true;
                while flag
                    temp = dir(fullfile(roi_path,'*_predictions.h5'));
                    temp2 = dir(fullfile(roi_path,'*.h5'));
                    if size(temp,1)>0 && size(temp2,1) == 2*size(temp,1)
                        flag = false;
                    else
                        fprintf('wating...');
                        pause(300);
                    end
                end
                H5_name = {temp.name};
                h5_num = size(H5_name, 2);
                if h5_num>0
                    [raw, labels, t_trace] = step2_instance_seg(roi_path, H5_name, h5_num, sample_num);
                    if ~is_test
                        trace_num = size(t_trace,1);
                        tic;
                        fprintf(['analysis\t',num2str(trace_num),'\t']);
                        % for n = 1:trace_num
                        %     step3_analysis(raw, labels, t_trace(n,:), view_path, sample_num, j);
                        % end
                        if trace_num == 1
                            step3_analysis(raw, labels, t_trace(1,:,:), view_path, j);
                        end
                        toc;
                    else
                        test_path = fullfile(view_path,sprintf('test%d', j-1));
                        if ~exist(test_path, 'dir')
                            mkdir(test_path);
                        end
                        write_test_data(labels, t_trace, test_path);
                    end
                    
                end
            end
        end
        if is_test
            test_montage(view_path);
        end
    end


end

