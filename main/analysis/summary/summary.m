function summary(data_path,max_num)
    temp = dir(fullfile(data_path,'v*'));
    View_name = {temp.name};
    view_num = size(View_name, 2);
    view_num = min(view_num,max_num);
    save_name = cell(4,1);
    save_name{1} = fullfile(data_path,'actin.csv');
    save_name{2} = fullfile(data_path,'claer.csv');
    save_name{3} = fullfile(data_path,'surface.csv');
    save_name{4} = fullfile(data_path,'tublin.csv');
    save_name{5} = fullfile(data_path,'kill.csv');
    for i = 1:view_num
        view_path = fullfile(data_path,View_name{i});
        fprintf([View_name{i},'\n'])
        
        write_summary(fullfile(view_path,'actin_result.csv'),save_name{1},19.6,true)
        write_summary(fullfile(view_path,'claer_result.csv'),save_name{2},1,true)
        write_summary(fullfile(view_path,'surface_result.csv'),save_name{3},0.0385,true)
        write_summary(fullfile(view_path,'tublin_result.csv'),save_name{4},1,false)
        write_summary(fullfile(view_path,'kill_result.csv'),save_name{5},1,true)
    end
end

function write_summary(file_name,save_name,factor,ismax)
    A = readmatrix(file_name)*factor;
    if ismax
        A(isnan(A)) = -inf;
        B = A(:,2:end);
        C = max(B,[],2);
    else
        A(isnan(A)) = inf;
        B = A(:,2:end);
        C = min(B,[],2);
    end
    fid = fopen(save_name, 'a+', 'n', 'utf8');
    str = '';
    for i = 1:size(C,1)
        if isinf(C(i))
            C(i)=0;
        end
%         if C(i)>0
%             C(i)=1;
%         end
        if i == 1
            str = [str,num2str(C(i))];
        else
            str = [str,',',num2str(C(i))];
        end
    end
    str = [str,'\n'];
    fprintf(fid, str);
    fclose(fid);
%     xlswrite(save_name,A,i);

end

