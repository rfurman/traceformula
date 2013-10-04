for N in range(2,100):
  print "--"
  print N
  for chi in DirichletGroup(N):
   if chi(-1)==1:
    sages = sum(vector(nf.coefficients(100)) for nf in CuspForms(chi,4).newforms(names="x"))
    mine = vector(TrTnew(n,N,4,chi) for n in range(1,101))
    diff = sages-mine
    if abs(diff)>0.001:
        print diff
        print chi
        print sages
        print mine
    else:
        print "Good"
