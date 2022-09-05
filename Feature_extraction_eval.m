% ----------------------------------------------------------------------- %
%    File_name: Calib.m
%    Programmer: Seungjae Yoo
%
%    Last Modified: 2022_09_03
%
% ----------------------------------------------------------------------- %
clear all
clc

data_labels = ['1','2', '3', '4', '5', '6','7', '8', '9'];
ms = [16,32,64];

for data_label = data_labels
    for m = ms


        FILENAME = strcat('/Users/goldenyoo/Library/Mobile Documents/com~apple~CloudDocs/BioCAS_prepare/BCIIV_2a_mat/A0',data_label,'E_mat');
        load(FILENAME);
        FILENAME = strcat('/Users/goldenyoo/Library/Mobile Documents/com~apple~CloudDocs/BioCAS_prepare/BCIIV_2a_mat/true_labels/A0',data_label,'E.mat');
        load(FILENAME);

        s(find(isnan(s))) = min(min(s));

        Class_1 = [];
        Class_2 = [];
        Class_3 = [];
        Class_4 = [];
        reject = [];

        i=1;
        event = 1;
        while i <=length(h.EVENT.TYP)
            if h.EVENT.TYP(i)== 783
                if classlabel(event) == 1
                    Class_1 = [Class_1 h.EVENT.POS(i)];
                elseif classlabel(event) == 2
                    Class_2 = [Class_2 h.EVENT.POS(i)];
                elseif classlabel(event) == 3
                    Class_3 = [Class_3 h.EVENT.POS(i)];
                else
                    Class_4 = [Class_4 h.EVENT.POS(i)];
                end

                event = event + 1;
                i = i + 1;
            else
                i = i + 1;
            end
        end

        for i = 1:length(reject)
            k = reject(1,i);
            cc1 = find(Class_1 == k+500);
            cc2 = find(Class_2 == k+500);
            cc3 = find(Class_3 == k+500);
            cc4 = find(Class_4 == k+500);


            Class_1(cc1) = [];
            Class_2(cc2) = [];
            Class_3(cc3) = [];
            Class_4(cc4) = [];
        end

        for j = length(Class_1):-1:1
            Class_1 = [Class_1 Class_1(j)+5  Class_1(j)+10 Class_1(j)+15 Class_1(j)+20 Class_1(j)+25];
        end

        for j = length(Class_2):-1:1
            Class_2 = [Class_2 Class_2(j)+5  Class_2(j)+10 Class_2(j)+15 Class_2(j)+20 Class_2(j)+25];
        end

        for j = length(Class_3):-1:1
            Class_3 = [Class_3 Class_3(j)+5  Class_3(j)+10 Class_3(j)+15];
        end

        for j = length(Class_4):-1:1
            Class_4 = [Class_4 Class_4(j)+5  Class_4(j)+10 Class_4(j)+15];
        end

        cue_offset = 313; % 1.25sec * 250Hz
        imagery_duration  = floor(2.75*250); % 3sec * 250Hz
        rear_offset = cue_offset + imagery_duration;


        channel_selection = [1:22];

        for i = 1:length(Class_1)
            group_1(i,:,:) = s(Class_1(i)+cue_offset:Class_1(i)+rear_offset,channel_selection);
            [K1(:,:,i), A1(:,:,i)] = my_SAX(squeeze(group_1(i,:,:)),m);
            Y1(i,1) = 1;
        end

        for i = 1:length(Class_2)
            group_2(i,:,:) = s(Class_2(i)+cue_offset:Class_2(i)+rear_offset,channel_selection);
            [K2(:,:,i), A2(:,:,i)] = my_SAX(squeeze(group_2(i,:,:)),m);
            Y2(i,1) = 2;
        end

        for i = 1:length(Class_3)
            group_3(i,:,:) = s(Class_3(i)+cue_offset:Class_3(i)+rear_offset,channel_selection);
            [K3(:,:,i), A3(:,:,i)] = my_SAX(squeeze(group_3(i,:,:)),m);
            Y3(i,1) = 3;
        end

        for i = 1:length(Class_4)
            group_4(i,:,:) = s(Class_4(i)+cue_offset:Class_4(i)+rear_offset,channel_selection);
            [K4(:,:,i), A4(:,:,i)] = my_SAX(squeeze(group_4(i,:,:)),m);
            Y4(i,1) = 4;
        end
        
        K1 = permute(K1,[2 1 3]);
        K2 = permute(K2,[2 1 3]);
        A1 = permute(A1,[2 1 3]);
        A2 = permute(A2,[2 1 3]);

        save_file_name = strcat('/Users/goldenyoo/Library/Mobile Documents/com~apple~CloudDocs/BioCAS_prepare/BCIIV_2a_mat/myData/Aug/Eval_data_',data_label ,'_chop_',int2str(m),'.mat')
        save(save_file_name,'K1','K2','A1','A2','Y1','Y2');
        clear group_1 group_2 group_3 group_4 K1 K2 K3 K4 A1 A2 A3 A4 Y1 Y2 Y3 Y4
    end
end
%%
function n_signal = my_normalization(s)
Mean = mean(s);
Std = std(s);

n_signal = (s - Mean)./Std;
end


function [K, A] = my_SAX(s,m)
n_signal = my_normalization(s);
q = floor(length(n_signal)/m);
K = [];
A = [];
for j = 1:m
    t = q*(j-1):q*j;
    v = n_signal(t+1,:);
    v_bar = mean(v);
    t_bar = mean(t);

    k_son = (t - t_bar)*(v - v_bar);
    k_mom = sum((t - t_bar).^2);

    k = k_son / k_mom;
    b = v_bar - k*t_bar;

    a = k*t_bar + b;

    K = [K k'];
    A = [A b'];
    %     A = [A a'];
end

end