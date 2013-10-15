#!/usr/bin/python
f = open('poly.txt','w');
numEvents = 100000;
sum = 0;
for i in range(100):
    sum += i*i*i*i*i;

for i in range(100):
    i5 = i*i*i*i*i;
    i5 = (i5*numEvents)/sum;
    round(i5);
    str = '{0} {1} '.format(i,i5);
    f.write(str+'\n');
f.close();
