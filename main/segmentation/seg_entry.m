clear
clc
ckpt = 'I:\CAR_T\code\segment\train_3-4\ckpt\best_checkpoint.pytorch';
data_path = 'I:\CAR_T\20230930_pd\20231017';
thread_num = 1;
% parpool(thread_num);
if thread_num > 1
    parfor i = 0:thread_num-1
        run_3dunet(data_path, ckpt, i, thread_num);
    end
else
    run_3dunet(data_path, ckpt, 0, thread_num);
end
