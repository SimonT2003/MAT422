% MAT 275 - Computer Lab 3
% By Simon Tran

%% Chapter 3 Questions

%% Problem #2 - Example 2
% This calculates the Wronskian of 3 functions
clear all
syms x
f1=exp(x);
f2=exp(-x);
f3=exp(2*x);
diff([f1 f2 f3], x) % 2nd row of Wronskian
A=[f1,f2,f3; diff([f1,f2,f3], x); diff([f1,f2,f3], x, 2)]; % 3x3 matrix
Wronsk=det(A);
subs(Wronsk, x, 0) % this calculates W(0)


%% Problem #3 - Example 3
% Complex Numbers
clear all
syms r1 r2 z x A1 A2

r1=2+3*1i;
real(r1)
imag(r1)
conj(r1)
r2=1/(4+5*1i);
exp(3+4*1i)
exp(pi*1i)
z=1-3*1i;
R=abs(z);
theta=angle(z);
R.*exp(1i*theta);
f(x)=x^4;
solve(f(x)-pi, x)
p(x)=x^2+6*x+25;
solve(p(x),x)
subs(p(x),x,1i)
subs(p(x),x,2+1i)
[solA1,solA2]=solve((2+1i)*A1==3, A1+3*A2*1i==0);


%% Problem #4 - Example 4
clear all
x0=1; xf=5;
y0=[2,0];
[x1,y1]=RK4(@Ch3NumExample1, [x0,xf],y0,0.05);
subplot(2,1,1), plot(x1,y1(:,1))
xlabel('x'); ylabel('y')
subplot(2,1,2), plot(x1,y1(:,2))
xlabel('x'); ylabel('y')
[x1(end), y1(end,:)] % this shows the last entry of vector x1 and 
% matrix y1. Note that we set xf=5.
figure % calls a new blank figure window
x0=1; xf=5;
y0=[1,0,0,1];
[x2,y2]=RK4(@Ch3NumExample2,[x0,xf],y0,.05);
plot(x2,y2)
legend('y(1)','y(2)','y(3)','y(4)')
xlabel('x');ylabel('y-values')
[x2(end),y2(end,:)]
figure
[x2,y2]=ode45(@Ch3NumExample2, [x0,xf],y0); % or use ode45
plot(x2,y2)
legend('y(1)','y(2)','y(3)','y(4)')
xlabel('x');ylabel('y-values')
[x2(end),y2(end,:)]


%% Problem #5 - Example 5
clear all
% Next three lines do NOT require the Symbolic Math Toolbox
p=[1 0 3 0 1]; % these are the coefficients in (v)
% if a coefficient is 0, you must put 0 in its position
roots(p)
% Now we use the Symbolic Math Toolbox:
syms x r
f(r)=r^2+3*r-4;
solve(f(r), r)
g(r)=r^4-4*r^3+6*r^2-4*r-15;
solve(g(r),r)
h(r)=r^4+r^3-2*r-1;
solve(h(r),r); % you can also do solve(h(r),r,'MaxDegree',4) to see the full answer
F(r)=r^3+r-3;
Fsoln=solve(F(r),r);
Fsoln(1) % This is the first root
double(Fsoln(1))
double(Fsoln)
G(r)=r^5+3*r^2-1;
Gsoln=solve(G(r),r);


%% Problem 9
clear all
syms x

% Part (a)
f1=exp(x);
f2=x*exp(x);
f3=x^2;
diff([f1 f2 f3], x) % 2nd row of Wronskian
A=[f1,f2,f3; diff([f1,f2,f3], x); diff([f1,f2,f3], x, 2)]; % this is a 3x3 matrix b/x we are given three functions
Wronsk=det(A);
subs(Wronsk, x, 0) % this calculates W(0)

% Part (b)
f4=sin(4*x)+cos(4*x);
f5=cos(4*x)-sin(4*x);
diff([f4 f5], x) % 2nd row of Wronskian
B=[f4,f5; diff([f4,f5], x)]; % this is 2x2 matrix b/c we are given two functions
Wronsk1=det(B);
subs(Wronsk1, x, 0) % this calculates W(0)


%% Problem 10
clear all 
syms r1 theta1 r2 theta2 r3 theta3 r4 theta4

% Part (a)
z1=2-3*1i;
real(z1)
imag(z1)
r1 = abs(z1); % this finds the r value
theta1 = angle(z1); % this finds the theta value
r1.*exp(1i*theta1); % this is in the form re^iθ

% Part (b)
z2=2i;
real(z2)
imag(z2)
r2 = abs(z2);
theta2 = angle(z2);
r2.*exp(1i*theta2);

% Part (c)
z3=-1-5i;
real(z3)
imag(z3)
r3 = abs(z3);
theta3 = angle(z3);
r3.*exp(1i*theta3);

% Part (d)
z4=3+1i;
real(z4)
imag(z4)
r4 = abs(z4);
theta4 = angle(z4);
r4.*exp(1i*theta4);


%% Problem 12a & 12c
syms z
p(z)=z^2+4*z+8;

% 12a
subs(p(z), z, 2i) % z = 2i

% 12c
subs(p(z), z, -1+4i) % z = -1 + 4i


%% Problem 16
clear all
x0=0; xf=20;
y0=[1, 2]; % IVP: y(0)=1, y'(0)=2
[x1,y1]=RK4(@Ch3Num16,[x0,xf],y0,0.1); % h=0.1 step size
plot(x1,y1)
legend('y(1)','y(2)')
xlabel('x'); ylabel('y')
[x1(end),y1(end,:)]

%% Problem 20
clear all
% Let's use the Symbolic Math Toolbox
syms r

% Part (a)
f(r) = r^2+r+1;
solve(f(r),r)

% Part (b)
g(r) = r^3+8*r^2+37*r+50;
solve(g(r), r)

% Part (c)
h(r) = r^4-4*r^3-2*r^2+36*r-63;
solve(h(r),r)


%% End of Code