clear
clc
data_path = 'I:\CAR_T\20230812_hd\data';
file_name = cell(4,1);
file_name{1} = fullfile(data_path,'actin.csv');
file_name{2} = fullfile(data_path,'claer.csv');
file_name{3} = fullfile(data_path,'surface.csv');
file_name{4} = fullfile(data_path,'kill.csv');
file_name{5} = fullfile(data_path,'tublin.csv');

C = [];
for i = 1:5
    A = readmatrix(file_name{i});
    C = cat(2,C,A(:));
end
% D = [];
% for i = 1:size(C,1)
%     if C(i,5)~=0
%         D = cat(1,D,C(i,:));
%     end
% end
[R,P] = corrcoef(C,'Rows','complete');