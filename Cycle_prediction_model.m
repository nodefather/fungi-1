clc
clear
close all
rng(10)

%% 加载数据
data = xlsread('\Fungal_trait_data.csv','D2:H35');
extension = xlsread('\extensionrate_moisture_temp.xlsx');

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


%% 建立模型 
Yt = data_in(:,37);
X1 = data_in(:,1:36);
X = tanh(X1);%增加非线性
W_out = pinv(X'*X+(1e-8).*eye(size(X,2)))*(X'*Yt);  
% figure
% bar(W_out);
%% 利用竞争规则更新模型菌株数量状态
%设定竞争排位对生长速率影响因子为rank_extend
%两两竞争关系，组合可以得到整体关系
%超参数设置
competitive_diff = 0.8;
Environment_diff = 1;

for i = 1:34
   for j = 1:34
       rank_extend(i,j) = competitive_diff*(fungi{j}.ranking-fungi{i}.ranking);%第j号菌株对第i号菌株的影响
   end
end
%% 验证模型(循环部分)
%给定环境时间序列
%时间序列长度
length = 365;
% 给定温度使-30到50℃之间
for x = 1:length
    Temperature(x,1) = 25*sin(x/30)+20;
    %给定湿度是-10到0Mpa之间
    Moisture(x,1) = 2.5*cos(x/30)-3;
end
% % 查看温度湿度随机曲线
% % figure
% % x = 1:length;
% % y = 1:length;
% % yyaxis left
% % plot(x,Temperature,'r');
% % xlabel('Time(day)');
% % ylabel('Temperature(℃)');
% % yyaxis right
% % plot(x,Moisture,'b');
% % ylabel('Moisture(Mpa)');
% % title('Moisture and Temperature curves');
% figure
% x = 1:length;
% y = 1:length;
% plot(x,Temperature);
% xlabel('Time(day)');
% ylabel('Temperature(℃)');
% title('Temperature curve');
% figure
% plot(x,Moisture);
% xlabel('Time(day)');
% ylabel('Moisture(Mpa)');
% title('Moisture curve');
% for x = 1:length
%     Temperature(x,1) = 20;
%     给定湿度是-10到0Mpa之间
%     Moisture(x,1) = -2.5;
% end

%取一组菌株均存在，且初始量相等的数据
in(1,1:34) = 1;
total = [];
extend_total = [];
for k =1:length
%湿度从Moisture中顺序选取
in(1,35) = Temperature(k,1);
%温度从Temperature中顺序选取
in(1,36) = Moisture(k,1);

%输入到模型中
Y(k,1) = tanh(in) * W_out;%得到分解率
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
total = [total;in];
extend_total = [extend_total;fungi{i}.extend_real];
end
%% 评价模型


figure
x = 1:length;
for i = 1:34
    plot(x,total(:,i));
    hold on
end
xlabel('Time(day)');
ylabel('mass(μg)');
title('The mass curve of fungal isolates(Environment and competition)');
