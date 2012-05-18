import csv
import matplotlib.pyplot as plt

time_point = 0
text = 'psbA time point ' + str(time_point+1)
gene = [413873,414955]
tss = 413826-gene[0]
#text = 'purF time point ' + str(time_point+1)
#gene = [4596, 6077]

walk = 1000

spread = [gene[0] - walk, gene[1] + walk]
px = []
nx = []
py = []
ny = []

f = open("GSE29264-GPL13535_series_matrix.txt")
while f.readline()[0]=='!':
  pass
val = list()
rows = csv.reader(f, delimiter='\t')
for r in rows:
  if r[0][0] == '!':
    break
  val.append([0]*len(r))
  for i in range(1,len(r)):
    val[-1][i-1] = float(r[i])
print 'Microarray imported'

f = open("GPL13535-11452.txt")
while f.readline()[0]=='#':
  pass
rows = csv.reader(f, delimiter='\t')
data = []
for r in rows:
  if len(r) > 0:
   data.append(r)
print 'Mapping imported'

for i in range(0,len(data)):
  start = int(data[i][4])
  stop = int(data[i][5])
  strand = data[i][3]
  db = data[i][2]
  type = data[i][1]
  if (start > spread[0]) and (stop < spread[1]) and (strand == '+') and (db == 'CP000100') and (type == 'chromosome'):
    s = 0
    for j in range(-5, 1):
      s += val[int(data[i+j][0])-1][time_point]
    s /= 5
    #~ s = val[int(data[i][0])-1][time_point]
    for j in range(0,12):
      px.append(start + j - gene[0])
      py.append(s)
  if (stop > spread[0]) and (start < spread[1]) and (strand == '-') and (db == 'CP000100')  and (type == 'chromosome'):
    s = 0
    for j in range(-5, 1):
      s += val[int(data[i+j][0])-1][time_point]
    s /= 5
    #~ s = val[int(data[i][0])-1][time_point]
    for j in range(0,12):
      nx.append(stop + j - gene[0])
      ny.append(s)
f.close()
print 'Data generated'

plt.plot(px,py,'r', label='Plus')
plt.plot(nx,ny,'b', label='Minus')
plt.fill_between([0,gene[1]-gene[0]], -1, 1, color='yellow', alpha=0.3)
plt.fill_between(px, 0, py, color='r', alpha=0.5)
plt.fill_between(nx, 0, ny, color='b', alpha=0.5)    
plt.axvline(tss)
plt.axvline()
plt.axvline(gene[1]-gene[0])
plt.gca().axis([spread[0]-gene[0], spread[1]-gene[0], -1, 1])
plt.xlabel('Genome position')
plt.ylabel('Expression Level')
plt.legend()
plt.title(text)
#plt.show()
plt.savefig(text + '.svg')
