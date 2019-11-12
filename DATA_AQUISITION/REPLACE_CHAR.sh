##### Python CSV read
### this script opens up a list of snotel files and changes ,s to .s
############################################################

import csv

x = []
with open ("list.csv") as myfile:
    data = myfile.read().split('\n')

length = len(data)

for i in range(length): 
    print(data[i]) 
    name = data[i]
    text = open(name, "r")
    text = ''.join([i for i in text]) \
    .replace(",", ".")
    x = open("output.csv","w")
    x.writelines(text)
    x.close()
    
    
    