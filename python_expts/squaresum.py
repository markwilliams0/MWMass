n=int(input('Enter a positive integer: '))
sum=0

for i in range(n+1):
    sum=sum+i**2
    
print 'The sum of integers from 1 to {} squared is {}'.format(n,sum)