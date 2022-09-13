import numpy as np

#calculem pi aaa

n = 1000000
sum = 0
for i in range(n-1):
    sum = sum + 1/(i+1)**2

pi = np.sqrt(6*sum)
print(f"pi = {pi}")
