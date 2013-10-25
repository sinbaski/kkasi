clear all
clc

syms x b;
assume(b > 0);

r = int('1/x^(3/2) * 1/(x + b)^(3/2)', x, -Inf, Inf);

