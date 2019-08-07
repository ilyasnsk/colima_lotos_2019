function int_best(a)

i1=int(a)
i2=i1+1
i3=i1-1

a1=i1
a2=i2
a3=i3
int_best=i1
if(abs(a2-a).lt.abs(a1-a)) int_best=i2
if(abs(a3-a).lt.abs(a1-a)) int_best=i3

return
end