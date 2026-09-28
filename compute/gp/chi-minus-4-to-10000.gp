default(parisize, "4G");
default(realprecision, 38);
L = lfuncreate(-4);
z = lfunzeros(L, 10000, 32);
for(i=1,#z, print(z[i]));
quit;
