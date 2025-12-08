function fout=f(x,Bi)
fout=((x*besselj(1,x))/(besselj(0,x)))-Bi;
end