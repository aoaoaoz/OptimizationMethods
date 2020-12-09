function y = dqdrtic(x)
%RPR 此处显示有关此函数的摘要
%   此处显示详细说明
y = 0; c = 100; d = 100;
for i=1:length(x)-2
    y = y + x(i).^2 + x(i+1).^2 *c + x(i+2).^2*d;
end
end