default(parisize, "2G");
G = znstar(5, 1);
L = lfuncreate([G, 2]);
print("# q=5 conrey=2 order=", charorder(G, 2), " chi(2)=", chareval(G, 2, 2, [I, 4]));
z = lfunzeros(L, 3000, 16);
for(i=1, #z, print(z[i]));
quit;
