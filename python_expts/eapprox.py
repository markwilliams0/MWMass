def fact(x):
    if x==0:
        return 1
    else:
        return float(x*fact(x-1))

n=int(input('Enter a positive integer: '))
if n%100==11 or n%100==12 or n%100==13:
    str='th'
elif n%10==1:
    str='st'
elif n%10==2:
    str='nd'
elif n%10==3:
    str='rd'
else:
    str='th'

est=0
for i in range(n+1):
    est=est+(1/(fact(i)))
print 'The value of e, to a {}{} order approximation, is {}'.format(n,str,est)