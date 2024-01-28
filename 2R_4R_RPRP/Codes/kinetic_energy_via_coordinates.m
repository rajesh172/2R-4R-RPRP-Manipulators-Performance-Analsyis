syms t th1(t) th2(t) r1(t) r2(t) r1_d(t) r2_d(t) r1_dd(t) r2_dd(t) th1_d(t) th2_d(t) th1_dd(t) th2_dd(t) x4(t) y4(t)
syms m1 m2 m3 m4 l1 l2
x4(t) = (l1+r1)*cos(th1) +(l2+r2)*cos(th2);
y4(t) = (l1+r1)*sin(th1) +(l2+r2)*sin(th2);

L = 0.5*m4*(diff(x4(t),t))^2 + (diff(y4(t),t))^2;
new_L = subs(L,{diff(th1(t),t),diff(th2(t),t),diff(th1_d(t),t),diff(th2_d(t),t),diff(r1(t),t),diff(r2(t),t),diff(r1_d(t),t),diff(r2_d(t),t)},{th1_d,th2_d,th1_dd,th2_dd,r1_d,r2_d,r1_dd,r2_dd})

% This will provide the kinetic energy of the mass if we know it's
% co-ordinates.





















