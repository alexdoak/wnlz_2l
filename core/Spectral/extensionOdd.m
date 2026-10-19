function [u] = extensionOdd(u)
u = [u; -u(end-1:-1:2)];
