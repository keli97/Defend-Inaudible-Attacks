function [yn,W,en]=LMS(xn,dn,M,mu,itr)
% LMS(Least Mean Squre)�㷨
% �������:
%     xn   ������ź�����      (������)
%     dn   ����������Ӧ����    (������)
%     M    �˲����Ľ���        (����)    �˲����Ľ���������ָ����г���Ĵ���,�����Խ�ߣ��˲�Ч����Խ��
%     mu   ��������(����)      (����)    Ҫ�����0,С��xn����ؾ����������ֵ�ĵ���    
%     itr  ��������            (����)    Ĭ��Ϊxn�ĳ���,M<itr<length(xn)
% �������:
%     W    �˲�����Ȩֵ����     (����)
%          ��СΪM : itr,
%     en   �������(itr : 1)    (������)  
%     yn   ʵ���������         (������)

% ������������Ϊ4����5��
if nargin == 4                 % 4��ʱ�ݹ�����Ĵ���Ϊxn�ĳ��� 
    itr = length(xn);
elseif nargin == 5             % 5��ʱ����M<itr<length(xn)
    if itr>length(xn) || itr<M
        error('��������������С!');
    end
else
    error('������������ĸ���!');
end

% ��ʼ������
en = zeros(itr,1);             % �������,en(k)��ʾ��k�ε���ʱԤ�������ʵ����������
W  = zeros(M,itr);             % ÿһ�д���һ����Ȩ����,ÿһ�д���-�ε���,��ʼΪ0

% ��������
for k = M:itr                  % ��k�ε���
    x = xn(k:-1:k-M+1);        % �˲���M����ͷ������
    
%     d = dn(k);
%     
%     en = sum((dn(k).^2-x.^2).^2);
%     [~,deltat] = min(en);
%     an = d(:,deltat)*(0.1:0.1:10);
%     en = sum((an.^2-x.^2).^2);
%     [~,deltat] = min(en);
%     x = x-an(:,deltat);
    y = W(:,k-1).' * x;        % �˲��������
    en(k) = (dn(k) - y) ;      % ��k�ε�������� 
    % �˲���Ȩֵ����ĵ���ʽ
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end

% ������ʱ�˲������������
yn = inf * ones(size(xn));
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;
end
