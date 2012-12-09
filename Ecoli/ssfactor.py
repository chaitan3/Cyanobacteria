from sklearn import svm, tree, naive_bayes, cross_validation
from sklearn.metrics import classification_report

#number of previous bases
nbp = 30
#number of bases afterwards
nba = 20
#number of k-Fold
nkf = 3
#rna-seq cutoff
cut = log2(4)

def transform(seqs):
  trans = []
  for seq in seqs:
    rnaseq = rna[seq]-rna[seq-1]
    mm = argmax(rnaseq)
    s = seq[0]+mm
    e = seq[-1]+1
    me = mean(rna[s:e])
    m = -inf
    for i in seq[:-5]:
      s = 0
      for j in range(0, 6):
        s += hmm[j][pseq[i+j]]
      if s > m:
        mi = i-seq[0]
        m = s
    trans.append([me*max(rnaseq), min(rnaseq)*int(argmin(rnaseq) < mm), m, mi-mm])
  return trans

data = []
ndata = []
for gene in genes:
  start = gene[0]
  end = gene[1]
  mrna = mean(rna[start:end])
  length = end - start
  if nbp <= start < len(rna)-nba and mrna > cut:
    #~ plot(range(start-nbp,start+nba), rna[start-nbp:start+nba]-rna[start-nbp-1:start+nba-1])
    #~ show()
    #~ print transform([arange(start-nbp, start+nba)])[0]
    data.append(transform([arange(start-nbp, start+nba)])[0])
    c = random()
    if c > 0.5:
      ndata.append(transform([arange(end-nbp,end+nba)])[0])
    else:
      r = int(random()*len(rna))
      while mean(rna[r-nbp:r+nba]) > 1:
        r = int(random()*len(rna))
      ndata.append(transform([arange(r-nbp,r+nba)])[0])
    
ss = [1 for i in data]
ssn = [0 for i in ndata]
ss.extend(ssn)
data.extend(ndata)
print 'Read positive & generated negative data'

c = svm.SVC()
#c = tree.DecisionTreeClassifier()
#c = naive_bayes.GaussianNB()
data = array(data)
ss = array(ss)
kfold = cross_validation.StratifiedKFold(ss, k=nkf, indices=False)
for train, test in kfold:
  c.fit(data[train], ss[train])
  sspred = c.predict(data[test])
  print 'Validation\n', classification_report(ss[test], sspred)
  
f = open('predicted.csv')
lines = f.readlines()
f.close()
pred = []
for line in lines:
  vals = line.split(',')
  if vals[1] == 'forward':
    start = int(vals[2])-1
    pred.append(arange(start-nbp,start+nba))
    print transform([pred[-1]])[0]
ypred = c.predict(transform(pred))
print 'Test\n', classification_report([1 for i in pred], ypred)
