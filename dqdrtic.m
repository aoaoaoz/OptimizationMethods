function y = dqdrtic(x)
%RPR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
y = 0; c = 100; d = 100;
for i=1:length(x)-2
    y = y + x(i).^2 + x(i+1).^2 *c + x(i+2).^2*d;
end
end