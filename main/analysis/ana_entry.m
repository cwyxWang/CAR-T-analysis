clear
clc
data_path = 'I:\CAR_T\20230930_pd\20231017';
t1 = clock;
step1_setup_analysis(data_path, 6, 0, false);
t2 = clock;
s = etime(t2,t1);

m = floor(s/60);
s = s-m*60;
h = floor(m/60);
m = m-h*60;
fprintf(['分析总耗时',num2str(h),'h',num2str(m),'m',num2str(s),'s\n'])