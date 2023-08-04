syms t th1(t) th2(t) r1(t) r2(t) r1_d(t) r2_d(t) r1_dd(t) r2_dd(t) th1_d(t) th2_d(t) th1_dd(t) th2_dd(t)
syms m1 m2 m3 m4 l1 l2
L1 = 1/2*(m1*l1^2)*(th1_d)^2;
L2 = 1/2*m2*((l1+r1)^2)*(th1_d^2) + 1/2*m2*(r1_d^2);
% L3 = 1/2*m3*(r1_d^2) + m3*l2*r1_d*th2_d*sin(th1-th2) + 1/2*m3*(l2)^2*(th2_d^2) + 1/2*m3*((l1+r1)^2)*(th1_d^2) + m3*l2*(l1+r1)*th1_d*th2_d*cos(th1-th2);
% L4 = 1/2*m4*(r1_d^2+r2_d^2) + m4*(l2+r2)*r1_d*th2_d*sin(th1-th2) + 1/2*m4*(l2)^2*(th2_d^2) + 1/2*m4*((l1+r1)^2)*(th1_d^2) + m4*(l2+r2)*(l1+r1)*th1_d*th2_d*cos(th1-th2) + m4*r2*(th2_d^2)*(l2+ r2/2) + m4*r1_d*r2_d*cos(th1-th2) - m4*(l1+r1)*r2_d*th1_d*sin(th1-th2);
L3 = m3/2*((sin(th1(t))*r1_d(t) + cos(th1(t))*th1_d(t)*(l1 + r1(t)) + l2*cos(th2(t))*th2_d(t))^2 + (l2*sin(th2(t))*th2_d(t) - cos(th1(t))*r1_d(t) + sin(th1(t))*th1_d(t)*(l1 + r1(t)))^2);
L4 = m4/2*((sin(th1(t))*r1_d(t) + sin(th2(t))*r2_d(t) + cos(th1(t))*th1_d(t)*(l1 + r1(t)) + cos(th2(t))*th2_d(t)*(l2 + r2(t)))^2 + (cos(th1(t))*r1_d(t) + cos(th2(t))*r2_d(t) - sin(th1(t))*th1_d(t)*(l1 + r1(t)) - sin(th2(t))*th2_d(t)*(l2 + r2(t)))^2);
L = L1 + L2 + L3 + L4;

tau1 = diff(diff(L,th1_d),t) - diff(L,th1);
F1 = diff(diff(L,r1_d),t) - diff(L,r1);
tau2 = diff(diff(L,th2_d),t) - diff(L,th2);
F2 = diff(diff(L,r2_d),t) - diff(L,r2);

new_tau1 = subs(tau1,{diff(th1(t),t),diff(th2(t),t),diff(th1_d(t),t),diff(th2_d(t),t),diff(r1(t),t),diff(r2(t),t),diff(r1_d(t),t),diff(r2_d(t),t)},{th1_d,th2_d,0,0,r1_d,r2_d,0,0})
new_F1 = subs(F1,{diff(th1(t),t),diff(th2(t),t),diff(th1_d(t),t),diff(th2_d(t),t),diff(r1(t),t),diff(r2(t),t),diff(r1_d(t),t),diff(r2_d(t),t)},{th1_d,th2_d,0,0,r1_d,r2_d,0,0});
new_tau2 = subs(tau2,{diff(th1(t),t),diff(th2(t),t),diff(th1_d(t),t),diff(th2_d(t),t),diff(r1(t),t),diff(r2(t),t),diff(r1_d(t),t),diff(r2_d(t),t)},{th1_d,th2_d,0,0,r1_d,r2_d,0,0});
new_F2 = subs(F2,{diff(th1(t),t),diff(th2(t),t),diff(th1_d(t),t),diff(th2_d(t),t),diff(r1(t),t),diff(r2(t),t),diff(r1_d(t),t),diff(r2_d(t),t)},{th1_d,th2_d,0,0,r1_d,r2_d,0,0});
 % we have put the value of all the double derivatives to be zero hence
 % this will provide the (Dynamics-Inertia) terms.

