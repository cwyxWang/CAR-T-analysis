data_path = 'I:\CAR_T\20230812_hd\data';
temp = dir(fullfile(data_path,'v*'));
View_name = {temp.name};
view_num = size(View_name, 2);
save_name = fullfile(data_path,'tublin.csv');
for i = 1:view_num
    view_path = fullfile(data_path,View_name{i});
    fprintf([View_name{i},'\n'])
    write_summary(fullfile(view_path,'tublin_result.csv'),save_name,1)

end

function write_summary(file_name,save_name,factor)
A = readmatrix(file_name)*factor;
A = A(:,2:end);
C = A(~isnan(A));

fid = fopen(save_name, 'a+', 'n', 'utf8');
str = '';
for i = 1:size(C,1)

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

