clc
clear
close all
rng(10)

%% 加载数据
data = xlsread('\Fungal_trait_data.csv','D2:H35');
extension = xlsread('\extensionrate_moisture_temp.xlsx');
T_cur = xlsread('Fungi_temperature_curves.csv');
Moi_cur = xlsread('Fungi_moisture_curves.csv');

%查看温湿度曲线
% figure
% for i = 1:35
%     plot(T_cur(5501*i-5500:5501*i,1),T_cur(5501*i-5500:5501*i,2));
%     hold on
% end
% xlabel('Temperature(℃)');
% ylabel('Extension Rate(mm/day)');
% title('The extension rate of some fungal isolates affected by temperature');

% figure
% for i = 1:35
%     plot(Moi_cur(501*i-500:501*i,1),Moi_cur(501*i-500:501*i,2));
%     hold on
% end
% xlabel('Moisture(Mpa)');
% ylabel('Extension Rate(mm/day)');
% title('The extension rate of some fungal isolates affected by moisture');
%% 处理数据
[a,b] = size(data);
%设定初始质量
%变量名为质量、湿度、温度
for i = 1:a
    data_in(i,i) = 0.5; %均设定为0.5
    data_in(i,a+1:a+b-2) = data(i,1:b-2);
end

%利用回归计算获取温湿度对每个菌种生长速率影响二维函数
%近似为非线性回归，拟合更加精确
%fungi{i}是含有真菌特性的元胞
for i = 1:a
    fungi{i}.X1 = extension(i:i+12,1:2);
    %归一化处理，归一化到[0,1]之间
    nMax = max(fungi{i}.X1);
    nMin = min(fungi{i}.X1);
    for j = 1:2
        fungi{i}.X1(:,j) = (fungi{i}.X1(:,j)- nMin(:,j))/(nMax(:,j)-nMin(:,j));
    end
    fungi{i}.X = ((-1).*fungi{i}.X1.^2+fungi{i}.X1);%假设该函数为二次函数非线性函数
    fungi{i}.Yt = extension(i:i+12,3);
    fungi{i}.extendWout = pinv((fungi{i}.X)'*(fungi{i}.X)+(1e-8).*eye(size(fungi{i}.X,2)))*((fungi{i}.X)'*(fungi{i}.Yt));%真菌生长速率随温湿度变化
    fungi{i}.density = data(i,b-1);%真菌密度
    fungi{i}.ranking = data(i,b);%真菌竞争排位
end

%绘制温湿度图像，湿度对应-5-0，温度对应10-40℃
for i = 1:10
    for j = 1:10
        ext_tem_moi(i,j) = ((-1).*[i/10,j/10].^2+[i/10,j/10]) * fungi{1}.extendWout;
%         ext_tem_moi(i,j) = [i/2-5,j*3+10] * fungi{1}.extendWout;
    end
end

% figure
% x = 0.1:0.1:1;
% y = 0.1:0.1:1;
% surfl(x,y,ext_tem_moi)
% xlabel('Moisture(Mpa)');
% ylabel('Temperature(℃)');
% zlabel('Extension Rate(mm/day)');
% title('Extension Rate affected by temperature and moisture');

%验证
% X = [-2,22];
% Y = tanh(X) * fungi{1}.extendWout;
%则得到Y是该温度和湿度下1号菌株的生长速率

%% 建立模型 
Yt = data_in(:,37);
X1 = data_in(:,1:36);
% %归一化
% nMax = max(X1);
% nMin = min(X1);
% for j = 1:36
%      X1(:,j) = (X1(:,j)- nMin(:,j))/(nMax(:,j)-nMin(:,j));
% end
X = (X1).^2;%增加非线性
W_out = pinv(X'*X+(1e-8).*eye(size(X,2)))*(X'*Yt);  
figure
bar(W_out);
xlabel('Features');
ylabel('Influence ability');
title('Effect of features on wood decomposition rate');

%% 利用竞争规则更新模型菌株数量状态
%设定竞争排位对生长速率影响因子为rank_extend
%两两竞争关系，组合可以得到整体关系
%超参数设置
competitive_diff = 0;
Environment_diff = 0.5;

for i = 1:34
   for j = 1:34
       rank_extend(i,j) = competitive_diff*(fungi{j}.ranking-fungi{i}.ranking);%第j号菌株对第i号菌株的影响
   end
end


%该因子乘以菌株质量即可得到对其他菌株生长速率的负影响
for k = 1:a
    %对每次数据的计算方法
%     %归一化
    for j = 1:2
%         fungi{i}.X1(:,j) = (fungi{i}.X1(:,j)- nMin(:,j))/(nMax(:,j)-nMin(:,j));
        extend_in(:,j) = (data_in(k,37-j)-nMin(:,j))/(nMax(:,j)-nMin(:,j));
    end
    
    for i = 1:34
        for j = 1:34
            %第i个菌株的生长速率
            fungi{i}.extend_real =  Environment_diff*((-1).*extend_in.^2+extend_in)*fungi{i}.extendWout - rank_extend(i,j)* data_in(k,j);
        end
    end
    %质量计算
    for i = 1:34
        %第i个菌株的菌落半径
        r = sqrt(data_in(k,i)/fungi{i}.density/pi);
        %新增面积计算函数
        fun = @(x) 2*pi*x;
        %质量更新,初始质量为0则不生长
        if data_in(k,i)>0
            if r+fungi{i}.extend_real>0
                data_in(k,i)=data_in(k,i) + 0.001*fungi{i}.density*integral(fun,r,r+fungi{i}.extend_real);
            end
            if r+fungi{i}.extend_real<=0
                data_in(k,i)=0;
            end
        end
        if data_in(k,i)<=0
            data_in(k,i) = 0;
        end
    end
end


%% 验证模型
%取一组菌株均存在，且初始量相等的数据
in(1,1:34) = 0.5;
%湿度从-5-0之间随机选取
in(1,36) = -2.5;
%温度从10-40℃之间随机选取
in(1,35) = 25;

% for j = 1:36
%      in(:,j) = (in(:,j)- nMin(:,j))/(nMax(:,j)-nMin(:,j));
% end
%输入到模型中
Y = (in).^2 * W_out+55;%得到分解率

%更新质量参数
for j = 1:2
    extend_in(:,j) = (in(1,37-j)-nMin(:,j))/(nMax(:,j)-nMin(:,j));
end
    
for i = 1:34
    for j = 1:34
    %第i个菌株的生长速率
        fungi{i}.extend_real =  Environment_diff*((-1).*extend_in.^2+extend_in)*fungi{i}.extendWout - rank_extend(i,j)* in(1,j);
    end
end
    %质量计算
for i = 1:34
        %第i个菌株的菌落半径
        r = sqrt(in(1,i)/fungi{i}.density/pi);
        %新增面积计算函数
        fun = @(x) 2*pi*x;
        %质量更新,初始质量为0则不生长
    if in(1,i)>0
        if r+fungi{i}.extend_real>0
            in(1,i)=in(1,i) + 0.001*fungi{i}.density*integral(fun,r,r+fungi{i}.extend_real);
        end
        if r+fungi{i}.extend_real<=0
            in(1,i) = 0;
        end
    end
    if in(1,i)<=0
        in(1,i) = 0;
    end
end

%% 评价模型
