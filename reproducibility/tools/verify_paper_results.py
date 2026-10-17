"""Read-only exact verification of the paper's core, using Python stdlib only.

No optimizer, network, personal config, session state or historical imports.
Default: verify frozen fixture and print a compact result. --emit PATH also
writes display data for MATLAB; that explicit mode is not read-only.
Proofs about function classes, smoothing and continuity remain analytic.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt, comb
import argparse, hashlib, json, platform, time

if not __debug__:
    raise RuntimeError('Run without -O: exact certificate assertions must remain enabled.')

ROOT=Path(__file__).resolve().parents[1]
ZERO=Q(0)
def decimal_bound(x, upper=False, digits=10):
    scaled=x*10**digits
    n=-((-scaled.numerator)//scaled.denominator) if upper else scaled.numerator//scaled.denominator
    sign='-' if n<0 else ''
    whole,part=divmod(abs(n),10**digits)
    result=f'{sign}{whole}.{part:0{digits}d}'
    assert (Q(result)>=x) if upper else (Q(result)<=x)
    return result

def qs(a): return [qs(x) for x in a] if isinstance(a,list) else Q(a)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),ZERO)
def mv(a,x): return [dot(row,x) for row in a]
def trans(a): return list(map(list,zip(*a)))
def radical(x):
    assert x>=0
    scale=10**60; n=isqrt(x.numerator*scale*scale//x.denominator)
    lo=Q(n,scale); hi=lo if lo*lo==x else Q(n+1,scale)
    assert lo*lo<=x<=hi*hi
    return lo,hi
def integeral(coeff): return sum((a/Q(i+1) for i,a in enumerate(coeff)),ZERO)
def absolute_integral(coeff,sub=64):
    """Bernstein convex-hull enclosure on [0,1], exact rational arithmetic."""
    degree=len(coeff)-1; upper=lower=ZERO
    for k in range(sub):
        x=Q(k,sub); h=Q(1,sub)
        power=[h**i*sum((coeff[j]*comb(j,i)*x**(j-i) for j in range(i,degree+1)),ZERO) for i in range(degree+1)]
        bern=[sum((power[j]*Q(comb(i,j),comb(degree,j)) for j in range(i+1)),ZERO) for i in range(degree+1)]
        lower+=h*abs(sum(bern,ZERO))/Q(degree+1)
        upper+=h*sum(map(abs,bern),ZERO)/Q(degree+1)
    assert lower<=upper
    return lower,upper
def cap_upper(p,s2,B,T):
    """Outward upper of coordinate support of ball intersect half-space."""
    level=sorted([p*x for x in T]); sl=radical(s2)[0]
    if B[0]>=level[1]: return T[1]
    if not(B[1]<level[0] and B[0]>-sl*T[0]): return T[1]
    b2min=ZERO if B[0]<=0<=B[1] else min(x*x for x in B)
    b2max=max(x*x for x in B)
    assert s2*T[0]**2-b2max>0
    return (max(p*x for x in B)+radical(s2-p*p)[1]*radical(s2*T[1]**2-b2min)[1])/s2

def verify(path):
    start=time.monotonic(); d=json.loads(path.read_text(encoding='utf-8-sig'))
    assert d['schema']==1 and d['wavelengthNm']==[380,780] and d['stepNm']==1
    C,F=qs(d['cmfsExact']),qs(d['spdsExact']); assert len(C)==len(F)==3
    n=len(C[0]); assert n==401 and all(len(x)==n for x in C+F)
    assert all(v>=0 for row in C+F for v in row)
    kernels=[]; T=[]; cs=[]
    for spd in F:
        block=[]
        for row in C:
            block.append([[row[i]*spd[i],(row[i+1]-row[i])*spd[i]+row[i]*(spd[i+1]-spd[i]),(row[i+1]-row[i])*(spd[i+1]-spd[i])] for i in range(n-1)])
        norm=100/sum((integeral(p) for p in block[1]),ZERO)
        block=[[[norm*x for x in p] for p in channel] for channel in block]
        kernels.extend(block); cs.append([sum((integeral(p) for p in channel),ZERO) for channel in block])
        for channel in block:
            row=[ZERO]*n
            for i,p in enumerate(channel):
                moment=sum((x/Q(k+2) for k,x in enumerate(p)),ZERO)
                row[i]+=integeral(p)-moment; row[i+1]+=moment
            T.append(row)
    assert all(c[1]==100 for c in cs)
    H=[[Q(1),-cs[0][0]/100,ZERO],[ZERO]*3,[ZERO]*3]
    hc=[mv(H,c) for c in cs]; h=sum(hc,[]); cn2=[dot(c,c) for c in cs]
    K2=max(dot(x,x)/c for x,c in zip(hc,cn2)); assert K2>0 and hc[0]==[ZERO]*3
    assert all(K2*dot(cs[i],cs[j])-dot(hc[i],hc[j])>=0 for i in range(3) for j in range(3))
    minor=[[row[i] for i in d['rankColumns']] for row in T]; iv=qs(d['rankInverseExact'])
    assert [[dot(row,col) for col in trans(iv)] for row in minor]==[[Q(i==j) for j in range(9)] for i in range(9)]
    l,u=qs(d['envelopeExact']); rho=Q(d['rhoExact']); betas={}; physical=None
    assert (l,u,rho)==(Q(-47,5000),Q(61,2500),Q(1,10))
    assert sorted(w['name'] for w in d['witnesses'])==['S0','S1','U0','U1']
    for w in d['witnesses']:
        z=qs(w['zExact']); beta=Q(w['betaExact']); r=Q(w['rhoExact'])
        assert w['slope']==w['name'].startswith('S')
        assert r==(rho if w['name'].endswith('1') else 0)
        assert len(z)==n and beta>0 and max(map(abs,z))<=Q(1,2)
        if w['slope']: assert all(l<=b-a<=u for a,b in zip(z,z[1:]))
        y=mv(T,z); assert y[:3]==[ZERO]*3
        e=[[a-beta*b for a,b in zip(y[3*j:3*j+3],hc[j])] for j in range(3)]
        assert all(r*r*beta*beta*K2*dot(cs[i],cs[j])-dot(e[i],e[j])>=0 for i in range(3) for j in range(3))
        betas[w['name']]=beta
        if w['name']=='S1':
            physical=[dict(illuminant=['D65','A','FL11'][j],target=list(map(float,[beta*a for a in hc[j]])),actualChange=list(map(float,y[3*j:3*j+3])),error=list(map(float,e[j])),absoluteErrorUpper=float(radical(dot(e[j],e[j]))[1]),globalRelativeErrorUpper=float(radical(dot(e[j],e[j])/(beta*beta*K2*cn2[j]))[1])) for j in range(3)]
    print('Exact inputs, rank, four witnesses and all-mixture feasibility passed.',flush=True)
    def setup(rec):
        mu=qs(rec['muExact']); flat=sum(mu,[])
        g=[[sum((flat[j]*kernels[j][i][k] for j in range(9)),ZERO) for k in range(3)] for i in range(n-1)]
        P=dot(flat,h); W=sum((radical(K2*cn2[j]*dot(mu[j],mu[j]))[1] for j in [1,2]),ZERO)
        assert P>0
        return mu,g,P,W
    mu,g,P,W=setup(d['supportU0'])
    U0upper=sum((absolute_integral(p)[1] for p in g),ZERO)/(2*P)
    assert U0upper==Q(d['supportU0']['upperExact']) and U0upper>=betas['U0']
    mu,g,P,W=setup(d['supportU1']); L=sum((absolute_integral(p)[1] for p in g),ZERO)/2
    assert L==Q(d['supportU1']['LExact']) and P-rho*W>0
    U1upper=L/(P-rho*W); uv=(L,P,W)
    print('Full measurable-class support upper bounds passed.',flush=True)
    mu,g,P,W=setup(d['supportS0']); eta=qs(d['supportS0']['etaExact'])
    dt=[-eta[0]]+[eta[i-1]-eta[i] for i in range(1,n-1)]+[eta[-1]]
    residue=[a-b for a,b in zip(mv(trans(T),sum(mu,[])),dt)]
    finite=(sum((u*a if a>=0 else l*a for a in eta),ZERO)+sum(map(abs,residue),ZERO)/2)
    peano=[[a/2+b/6+c/12,-a,-b/2,-c/3] for a,b,c in g]
    assert all(integeral(p)==0 for p in peano)
    correction=(u-l)*sum((absolute_integral(p)[1] for p in peano),ZERO)/2
    B0=(finite+correction)/P; assert B0==Q(d['supportS0']['upperExact'])
    smu,sP,sL=mu,P,P*B0
    mu,g,P,W=setup(d['supportS1']); pieces=qs(d['supportS1']['piecesExact'])
    assert len(pieces)==800 and pieces[0][0]==0 and sum(pieces[-1],ZERO)==0
    assert all(sum(a,ZERO)==b[0] for a,b in zip(pieces,pieces[1:]))
    residual=ZERO; signed=ZERO; absup=ZERO
    for j,p in enumerate(pieces):
        a,b,c=g[j//2]; x=Q(j%2,2)
        gp=[a+b*x+c*x*x,(b+2*c*x)/2,c/4]
        difference=[2*p[1]-gp[0],4*p[2]-gp[1],6*p[3]-gp[2]]
        assert difference[1:]==[ZERO]*2
        residual+=abs(difference[0])/4
        signed+=integeral(p)/2
        absup+=absolute_integral(p,d['supportS1']['subdivisions'])[1]/2
    assert P-rho*W>0
    S1upper=((u-l)*absup/2-(u+l)*signed/2+residual)/(P-rho*W)
    assert S1upper==Q(d['supportS1']['upperExact'])
    print('Full AC bounds: finite support + Peano, and cubic integration-by-parts passed.',flush=True)
    absolute=d['absolute']; vals={k:Q(v) for k,v in absolute['valuesExact'].items()}
    b0,b1=betas['S0'],betas['S1']; U0=betas['U0']; L,P,W=uv
    assert [b0,b1,B0,U0,L,P,W]==[vals[k] for k in ['envelopeBaselineWitness','envelopeTolerantWitness','envelopeBaselineUpper','unrestrictedBaselineLower','supportUpper','pairing','penaltyUpper']]
    k=(b1-b0)/(rho*b1); assert k==vals['blendKappa']
    ds=[b1-B0/(1-rho),S1upper-b0/(1-rho)]
    du=[betas['U1']-U0upper/(1-rho),U1upper-U0/(1-rho)]
    assert ds==qs(absolute['restrictedBonusEnclosureExact']) and du==qs(absolute['unrestrictedBonusEnclosureExact'])
    assert du[0]>0 and ds[0]>Q(3,2)*du[1] and B0<U0 and S1upper<betas['U1']
    gammaS=[b1-B0,S1upper-b0];gammaU=[betas['U1']-U0upper,U1upper-U0]
    assert gammaS[0]-gammaU[1]>Q('0.0355394472')
    # Correct relative-penalty proof: upper tolerant gap / lower exact gap.
    penalty_lower=1-(U1upper-b1)/(U0-B0)
    printed_d0=Q('0.3817259867')-Q('0.2171923873')
    printed_d1=Q('0.5317220254')-Q('0.4027278731')
    displayed_penalty=1-printed_d1/printed_d0
    assert U0-B0>=printed_d0>0 and 0<=U1upper-b1<=printed_d1
    assert displayed_penalty==Q(355394471,1645335994)>Q('0.2160')
    assert penalty_lower>Q('0.21600115')
    capacity_exact=[[U0,U0upper],[b0,B0],[betas['U1'],U1upper],[b1,S1upper]]
    printed_capacity=[[decimal_bound(lo),decimal_bound(hi,True)] for lo,hi in capacity_exact]
    assert printed_capacity==[['0.3817259867','0.3835151503'],['0.2101706137','0.2171923873'],['0.5127794981','0.5317220254'],['0.4027278731','0.4863981307']]
    printed_bonus=[[decimal_bound(lo),decimal_bound(hi,True)] for lo,hi in [du,ds]]
    assert printed_bonus==[['0.0866515534','0.1075820401'],['0.1614029983','0.2528752265']]
    assert ds[0]/du[1]>Q('1.50027828')
    D=U0-B0; q=qs(absolute['lowerComparisonPolynomialExact'])
    assert q==[b0*P-L+D*P,-b0*(P+W)+L*(1+k)-D*(W+k*P),b0*W-L*k+D*k*W]
    linear=[((1-rho)*b1+D)*P-L,L-((1-rho)*b1+D)*W]
    left,right=qs(absolute['closedOrderingIntervalExact'])
    assert (left,right)==(Q('0.056736'),Q('0.137809'))
    assert left<rho<right<1 and q[0]<0<q[1] and q[2]>0
    assert q[0]+q[1]*left+q[2]*left*left>0 and linear[1]<0 and linear[0]+linear[1]*right>0
    assert 1-k*rho>0 and P-W*right>0
    # Check separate R133 scope algebra, not new joint direction/norm claims.
    sc=d['scope']; HZ=[[ZERO]*3,[ZERO]*3,[ZERO,-cs[0][2]/100,Q(1)]]
    hz=[mv(HZ,c) for c in cs]; KZ2=max(dot(x,x)/c for x,c in zip(hz,cn2))
    assert KZ2==Q(d['zDirection']['K2Exact'])
    assert all(KZ2*dot(cs[i],cs[j])-dot(hz[i],hz[j])>=0 for i in range(3) for j in range(3))
    zw=d['zWitness']; zz=qs(zw['baselineZExact']); bz=Q(zw['baselineBetaExact'])
    assert max(map(abs,zz))<=Q(1,2) and all(l<=b-a<=u for a,b in zip(zz,zz[1:]))
    assert mv(T,zz)==[bz*x for x in sum(hz,[])]
    zu=d['zUpper']; zm=qs(zu['physicalMultipliersExact']); ze=qs(zu['derivativeMultipliersExact'])
    zg=[[sum((sum(zm,[])[j]*kernels[j][i][k] for j in range(9)),ZERO) for k in range(3)] for i in range(n-1)]
    zp=dot(sum(zm,[]),sum(hz,[]));zwp=sum((radical(KZ2*cn2[j]*dot(zm[j],zm[j]))[1] for j in [1,2]),ZERO)
    assert zp-rho*zwp>0
    zdt=[-ze[0]]+[ze[i-1]-ze[i] for i in range(1,n-1)]+[ze[-1]]
    zr=[a-b for a,b in zip(mv(trans(T),sum(zm,[])),zdt)]
    zfinite=sum((u*a if a>=0 else l*a for a in ze),ZERO)+sum(map(abs,zr),ZERO)/2
    zpeano=[[a/2+b/6+c/12,-a,-b/2,-c/3] for a,b,c in zg]
    assert all(integeral(p)==0 for p in zpeano)
    zcorrection=(u-l)*sum((absolute_integral(p)[1] for p in zpeano),ZERO)/2
    zupper=(zfinite+zcorrection)/(zp-rho*zwp)
    assert zupper==Q(zu['fullACBetaUpperExact'])
    assert zupper-bz/(1-rho)==Q(sc['zBetaBonusUpperExact'])
    for Ksq,key in [(K2,'xKEnclosureExact'),(KZ2,'zKEnclosureExact')]:
        kl,ku=qs(sc[key]);assert kl*kl<=Ksq<=ku*ku
    ax=Q(sc['xKEnclosureExact'][0])*ds[0]
    az=Q(sc['zKEnclosureExact'][1])*Q(sc['zBetaBonusUpperExact'])
    assert ax>Q('45.278')*az>0
    nq=[L*k*rho*(1-rho)-D*k*rho*P,b0*P*(1-rho)-L*(1-rho)+D*(P+k*rho*rho*W),-rho*W*(b0*(1-rho)+D)]
    assert nq==qs(sc['normComparisonQuadraticExact'])
    c=Q(sc['normClosedConditionUpperExact']);t=c*c
    assert c==Q('1.052918')
    assert nq[2]<0 and nq[1]+2*nq[2]<0 and sum(nq)>0
    assert nq[0]+nq[1]*t+nq[2]*t*t>0 and 1-k*rho>0 and P-rho*W*t>0
    response=d['response']; assert Q(response['supportUpperExact'])==sL and Q(response['betaThresholdExact'])==b1
    for j,row in enumerate(response['lamps'],1):
        other=3-j;s2=dot(smu[j],smu[j]); rb=[rho*x for x in radical(K2*cn2[j])]
        ob=[rho*x for x in radical(K2*cn2[other]*dot(smu[other],smu[other]))]
        Ccap=[sP-ob[1],sP-ob[0]];B=[sL-b1*Ccap[1],sL-b1*Ccap[0]];TT=[b1*x for x in rb]
        deficit=-cap_upper(smu[j][1],s2,B,TT)
        actualX=b1*hc[j][0]-cap_upper(-smu[j][0],s2,B,TT)
        assert deficit==Q(row['YDeficitLowerExact']) and actualX==Q(row['actualXIncrementLowerExact'])
        assert Ccap[0]>rb[1]*radical(s2)[1] and hc[j][0]>rb[1]
        assert actualX/(B0*hc[j][0])-1==Q(row['actualXGainLowerExact'])>0
        pen=radical(s2-smu[j][1]**2)[1] if smu[j][1]>=0 else radical(s2)[1]
        assert Ccap[0]-rb[1]*pen>0
        assert sL/(Ccap[0]-rb[1]*pen)==Q(row['noYLossBetaUpperExact'])
    assert Q(response['lamps'][0]['YDeficitLowerExact'])>0
    assert Q(response['lamps'][0]['noYLossBetaUpperExact'])<b1
    for j,(xlower,xupper,gain) in enumerate([('5.44012708','3.21589624','0.691636'),('1.77273816','1.28444742','0.380156')]):
        row=response['lamps'][j]
        assert Q(row['actualXIncrementLowerExact'])>Q(xlower)
        assert Q(row['baselineXIncrementUpperExact'])<Q(xupper)
        assert Q(row['actualXGainLowerExact'])>Q(gain)
    assert Q(response['lamps'][0]['YDeficitLowerExact'])>Q('0.24083547')
    assert Q(response['lamps'][0]['noYLossBetaUpperExact'])<Q('0.2679937276830')<b1
    assert Q('0.4027278731')<b1 # Rounded figure label is not the exact threshold.
    assert ax>Q('0.0156453680148') and az<Q('0.0003455399445')
    assert rho*t<1
    assert -7 in d['sourceSolverFlags'] # Preserve unsuccessful solver provenance.
    print('Absolute comparison, closed tolerance interval and response-cost certificates passed.',flush=True)
    result=dict(passed=True,python=platform.python_version(),elapsedSeconds=time.monotonic()-start,fixtureSha256=hashlib.sha256(path.read_bytes()).hexdigest(),
      verification='exact rational core certificates; analytic propositions examined separately',readOnlyDefault=True,newOptimizerCalls=0,
      capacities=dict(U0=[float(U0),float(U0upper)],U1=[float(betas['U1']),float(U1upper)],S0=[float(b0),float(B0)],S1=[float(b1),float(S1upper)]),
      bonuses=dict(S=list(map(float,ds)),U=list(map(float,du))),gamma=dict(S=list(map(float,gammaS)),U=list(map(float,gammaU))),
      printedCapacityBounds=printed_capacity,printedBonusBounds=printed_bonus,
      penaltyReduction=dict(exactLower=str(penalty_lower),displayedLower=str(displayed_penalty),lowerApprox=float(penalty_lower)),
      ratioLower=float(ds[0]/du[1]),closedToleranceInterval=[float(left),float(right)],values=absolute['values'],
      physical=physical,responseBounds=response['lamps'],scope=dict(alphaXLower=float(ax),alphaZUpper=float(az),factorLower=float(ax/az),normConditionUpper=float(c)),
      limitations=['Not proof-assistant verification','No smooth endpoint attainment','No numerical input-uncertainty radius','No perceptual or manufacturing assertion'])
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--fixture',type=Path,default=ROOT/'results/research/paper-fixture.json');ap.add_argument('--emit',type=Path);args=ap.parse_args()
    out=verify(args.fixture)
    if args.emit: args.emit.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:out[k] for k in ['passed','python','elapsedSeconds','fixtureSha256','ratioLower','closedToleranceInterval']},indent=2))
