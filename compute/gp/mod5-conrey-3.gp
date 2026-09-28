default(parisize, "2G");
G = znstar(5, 1);
L = lfuncreate([G, 3]);
print("# q=5 conrey=3 order=", charorder(G, 3), " chi(2)=", chareval(G, 3, 2, [I, 4]));
z = lfunzeros(L, 3000, 16);
for(i=1, #z, print(z[i]));
quit;
