import csv
from matplotlib.pylab import plot, show
from numpy import arange
from scipy import interpolate
	
folder = 'E-GEOD-18902'
file_start = 'GSM'
#18902
file_num = 468463
#14225
#file_num = 356401
file_end = '_sample_table.txt'
#18902
ts = 24
te = 84
#14225
#ts = 4
#te = 48
ti = 4

#kaiA 18902
#gene = 637799647
#kaiA 14225
#gene = 'syc0334-m_at'

#purF 18902
#gene = 637798409
#purF 14225
#gene = 'syc1507-h_at'


#psbA1 18902
gene = 637798830
#psbA1 14225
#gene = 'syc1103-h_at'

#psbA2 18902
#gene = 637799310
#psbA2 14225
#gene = 'syc0167-m_s_at'

s = []
t = arange(ts,te+1,ti)
for i in range(0, len(t)):
	r = csv.reader(open(folder + '/' + file_start + str(file_num+i) + file_end),delimiter='\t')
	s.append(dict())
	for row in r:
		s[i][row[0]] = row[1]
y = []
for i in s:
	y.append(i[str(gene)])
	
plot(t,y)

#Interpolation
sp = interpolate.interp1d(t, y,kind='cubic')
tn = arange(ts, te, 0.1)
yn = sp(tn)
plot(tn, yn)


show()		
