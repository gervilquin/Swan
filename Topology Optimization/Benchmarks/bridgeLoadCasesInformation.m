% bridgeLoadCasesInformation %
% MESH 100x20 elements
% a = 39:20:219;
% b = 40:20:220;
% c = 259:20:439;
% d = 260:20:440;
% f = 480:20:660;
% e = f-1;
% h = 700:20:880;
% g2 = h-1;
% j = 920:20:1100;
% i2 = j-1;
% l = 1140:20:1320;
% k = l-1;
% n = 1360:20:1540;
% m = n-1;
% p = 1580:20:1760;
% o = p-1;
% r = 1800:20:1980;
% q = r-1;
% g = struct();
% for i = 1:length(a)
% g.zero(2*i -1)  =  a(i);
% g.zero(2*i) = b(i);
% g.one(2*i-1)  = c(i);
% g.one(2*i)  = d(i);
% g.two(2*i -1)  = e(i);
% g.two(2*i)  = f(i);
% g.three(2*i -1)  = g2(i);
% g.three(2*i)  = h(i);
% g.four(2*i -1)  = i2(i);
% g.four(2*i)  = j(i);
% g.five(2*i -1)  = k(i);
% g.five(2*i)  = l(i);
% g.six(2*i-1)  = m(i);
% g.six(2*i)  = n(i);
% g.seven(2*i-1)  = o(i);
% g.seven(2*i)  = p(i);
% g.eight(2*i-1)  = q(i);
% g.eight(2*i)  = r(i);
% end
g.zerov = [0 -10  0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.onev = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.twov = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.threev = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.fourv = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.fivev = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.sixv = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.sevenv = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.eightv = [0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10     0   -10];
g.zero = [ 79    80   119   120   159   160   199   200   239   240   279   280   319   320   359   360   399   400   439   440];
g.one  = [519   520   559   560   599   600   639   640   679   680   719   720   759   760   799   800   839   840   879   880];
g.two  = [ 959         960         999        1000        1039        1040        1079        1080        1119        1120 1159 ...
    1160        1199        1200        1239        1240        1279        1280        1319        1320];
g.three = [1399        1400        1439        1440        1479        1480        1519        1520        1559        1560 ...
    1599        1600        1639        1640        1679        1680        1719        1720        1759        1760];
g.four = [1839        1840        1879        1880        1919        1920        1959        1960        1999        2000 ...
 2039        2040        2079        2080        2119        2120        2159        2160        2199        2200];
g.five = [2279        2280        2319        2320        2359        2360        2399        2400        2439        2440 ...
    2479        2480        2519        2520        2559        2560        2599        2600        2639        2640];
g.six = [2719        2720        2759        2760        2799        2800        2839        2840        2879        2880 ...
    2919        2920        2959        2960        2999        3000        3039        3040        3079        3080];
g.seven = [3159        3160        3199        3200        3239        3240        3279        3280        3319        3320 ...
    3359        3360        3399        3400        3439        3440        3479        3480        3519        3520];
g.eight = [3599        3600        3639        3640        3679        3680        3719        3720        3759        3760 ...
    3799        3800        3839        3840        3879        3880        3919        3920        3959        3960];