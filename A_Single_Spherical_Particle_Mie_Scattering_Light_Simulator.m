tic

clear all

lambda = 550 * 10^-9; % 欲代入之入射光波長, 單位 <m>
k0 = 2 * pi / lambda; % wave munber
k = k0*1.5; % 定義參數k
r = 50 * 10^-9; % 粒徑大小, 單位 <m>
m = 2/1.5; % 粒子折射率
x = k * r; % 定義參數x
z = m * x; % 定義參數z
i = sqrt(-1); % 定義虛數i
N = round( x + 4.05*x^(1/3) + 2 ); % 定義無窮級數的近似項數
Z = 0.001; % 粒子的Z座標位置 

Er = cos(pi/4); % 定義入射光水平分量
El = sin(pi/4); % 定義入射光垂直分量

hold on
grid on

% 以下計算過程請參考論文
for th = 0:1:90 % th 為散射光與Z軸的夾角, 1 表一度為一個間隔
    
    theta = th * pi / 180; % 將 th 換算成徑度
    R = sec(theta) * Z; % 粒子離座標原點的距離

    PI(0 + 1) = 0;
    PI(1 + 1) = 1; 
    TAU(0 + 1) = 0;
    u = cos(theta);

    for p = 2:N+1
        s = u * PI(p);
        t = s - PI(p - 1);
        PI(p + 1) = s + t + t / (p - 1);
        TAU(p) = (p - 1) * t - PI(p - 1);
    end
 
    for t=1:N  
        T(t) = TAU(t+1);
        P(t) = PI(t+1);
    end
   
    for n = 1:N
        D(n) = -n/z + phi(n-1, z)/phi(n, z);
        a(n) =  ( (D(n)/m + n/x) * phi(n, x) - phi(n-1, x) ) / ( (D(n)/m + n/x) * kersi(n, x) - kersi(n-1, x) );
        b(n) = ( (m*D(n) + n/x) * phi(n, x) - phi(n-1, x) ) / ( (m*D(n) + n/x) * kersi(n, x) - kersi(n-1, x) );   
    end

    for u = 1:N 
        s1(u) = ( 2*u + 1 ) / ( u*(u+1) ) * ( a(u)*P(u) + b(u)*T(u) );
        s2(u) = ( 2*u + 1 ) / ( u*(u+1) ) * ( b(u)*P(u) + a(u)*T(u) );
    end
    
    S1 = sum(s1);
    S2 = sum(s2);
    Eout = (exp(-i*k*R + i*k*Z) / (i*k*R)) * [ S1, 0; 0, S2 ] * [ Er; El ];

    X(th+1, :) = Eout(1, 1); % 紀錄散射光的水平分量
    Y(th+1, :) = Eout(2, 1); % 紀錄散射光的垂直分量

    ang1 = angle( Eout(2, 1)/Eout(1, 1) ) * 180 / pi; % 計算 phase
    ANG1(th+1, :) = ang1; % 紀錄 phase

    plot(th, ang1, 'bo'); % 作散射角與phase的關係圖
end

xlabel('散射角 <deg.>')
ylabel('相位差 <deg.>')

toc