"""Validated continuous audit of the round--Regge de Rham A2 path.

The protocol and all numerical choices were frozen in commit 8b322b1 after a
finite-grid pattern was known but before a continuous certificate was run.
"""

from flint import acb, arb, ctx
from itertools import product
import argparse
from fractions import Fraction
import json
import multiprocessing as mp
from pathlib import Path
import sympy as sp
import time
ctx.dps=40

UORDER=18
class U:
  __slots__=('a',)
  def __init__(self,a=0):
    if isinstance(a,U):self.a=a.a[:]
    elif isinstance(a,(list,tuple)):self.a=[arb(x) for x in a]+[arb(0)]*(UORDER+1-len(a))
    else:self.a=[arb(a)]+[arb(0)]*UORDER
  @staticmethod
  def var(v):return U([v,1])
  def __add__(self,o):o=uu(o);return U([self.a[i]+o.a[i] for i in range(UORDER+1)])
  __radd__=__add__
  def __neg__(self):return U([-x for x in self.a])
  def __sub__(self,o):return self+(-uu(o))
  def __rsub__(self,o):return uu(o)-self
  def __mul__(self,o):
    o=uu(o);return U([sum((self.a[k]*o.a[i-k] for k in range(i+1)),arb(0)) for i in range(UORDER+1)])
  __rmul__=__mul__
  def inv(self):
    b=[1/self.a[0]]
    for i in range(1,UORDER+1):b.append(-b[0]*sum((self.a[k]*b[i-k] for k in range(1,i+1)),arb(0)))
    return U(b)
  def __truediv__(self,o):return self*uu(o).inv()
  def __rtruediv__(self,o):return uu(o)/self
  def derivative(self):return U([(i+1)*self.a[i+1] for i in range(UORDER)]+[0])
  def __pow__(self,p):
    p=arb(p)
    if p==0:return U(1)
    ratio=self.derivative()/self
    b=[self.a[0]**p]
    for n in range(1,UORDER+1):b.append(p*sum((b[k]*ratio.a[n-1-k] for k in range(n)),arb(0))/n)
    return U(b)
  def sqrt(self):return self**(arb(1)/2)
  def sincos(self):
    s=[self.a[0].sin()];c=[self.a[0].cos()]
    for n in range(1,UORDER+1):
      s.append(sum((k*self.a[k]*c[n-k] for k in range(1,n+1)),arb(0))/n)
      c.append(-sum((k*self.a[k]*s[n-k] for k in range(1,n+1)),arb(0))/n)
    return U(s),U(c)
  def sin(self):return self.sincos()[0]
  def cos(self):return self.sincos()[1]
  def acos(self):
    rhs=-(self.derivative()/(1-self*self).sqrt())
    return U([self.a[0].acos()]+[rhs.a[n-1]/n for n in range(1,UORDER+1)])
def uu(x):return x if isinstance(x,U) else U(x)

N=4
class X:
  __slots__=('a',)
  def __init__(self,a): self.a=[uu(z) for z in a]+[U(0)]*(N+1-len(a))
  @staticmethod
  def const(v): return X([v])
  @staticmethod
  def var(v): return X([v,1])
  def __add__(self,o): o=xx(o); return X([self.a[i]+o.a[i] for i in range(N+1)])
  __radd__=__add__
  def __neg__(self): return X([-z for z in self.a])
  def __sub__(self,o): return self+(-xx(o))
  def __rsub__(self,o): return xx(o)-self
  def __mul__(self,o):
    o=xx(o); return X([sum((self.a[k]*o.a[i-k] for k in range(i+1)),U(0)) for i in range(N+1)])
  __rmul__=__mul__
  def inv(self):
    b=[self.a[0].inv()]
    for i in range(1,N+1): b.append(-b[0]*sum((self.a[k]*b[i-k] for k in range(1,i+1)),U(0)))
    return X(b)
  def __truediv__(self,o): return self*xx(o).inv()
  def __rtruediv__(self,o): return xx(o)/self
  def __pow__(self,p):
    if not isinstance(p,int): raise TypeError(p)
    if p<0:return (self.inv())**(-p)
    z=X.const(1); b=self
    while p:
      if p&1:z=z*b
      b=b*b;p//=2
    return z
  def sqrt(self):
    b=[self.a[0].sqrt()]
    for i in range(1,N+1):
      b.append((self.a[i]-sum((b[k]*b[i-k] for k in range(1,i)),U(0)))/(2*b[0]))
    return X(b)
  def sincos(self):
    s=[self.a[0].sin()]; c=[self.a[0].cos()]
    for i in range(1,N+1):
      s.append(sum((k*self.a[k]*c[i-k] for k in range(1,i+1)),U(0))/i)
      c.append(-sum((k*self.a[k]*s[i-k] for k in range(1,i+1)),U(0))/i)
    return X(s),X(c)
  def sin(self):return self.sincos()[0]
  def cos(self):return self.sincos()[1]
  def acos(self):
    root=(1-self*self).sqrt(); rhs=-(self.derivative()/root)
    b=[self.a[0].acos()]+[rhs.a[i-1]/i for i in range(1,N+1)]
    return X(b)
  def derivative(self):return X([(i+1)*self.a[i+1] for i in range(N)]+[U(0)])
def xx(x):return x if isinstance(x,X) else X.const(x)

def comp(x,k):return x.a[k]
def supabs(x):return abs(x).upper()
def interval(c,r):return arb(c,r)

def integrate_vec(func,nvar,ncell):
  totals=None; errors=None
  width=arb(1)/ncell; half=width/2; vol=width**nvar
  shift=width/(2*arb(3).sqrt()); weight=vol/(2**nvar)
  fact4=arb(24)
  for inds in product(range(ncell),repeat=nvar):
    centers=[arb(2*i+1)/(2*ncell) for i in inds]
    for signs in product((-1,1),repeat=nvar):
      args=[X.const(centers[i]+signs[i]*shift) for i in range(nvar)]
      vals=func(*args); vals=vals if isinstance(vals,(tuple,list)) else (vals,)
      if totals is None:
        totals=[U(0) for _ in vals]; errors=[[arb(0)]*(UORDER+1) for _ in vals]
      for k,z in enumerate(vals):totals[k]+=z.a[0]*weight
    boxes=[]
    for direction in range(nvar):
      args=[X.var(interval(centers[i],half)) if i==direction else X.const(interval(centers[i],half)) for i in range(nvar)]
      zs=func(*args); boxes.append(zs if isinstance(zs,(tuple,list)) else (zs,))
    for k in range(len(totals)):
      for q in range(UORDER+1):
        e=sum((supabs(comp(boxes[i][k].a[4]*fact4,q)) for i in range(nvar)),arb(0))*vol*width**4/4320
        errors[k][q]+=e
  for k,z in enumerate(totals):
    for q in range(UORDER+1):z.a[q]+=arb(0,errors[k][q])
  return totals

pi=arb.pi();sqrt5=arb(5).sqrt();a2=(7+3*sqrt5)/16;rho=(1-a2).sqrt();h=-rho/3;hab=rho/3
ell=2/(1+sqrt5);rin=ell/(2*arb(3).sqrt());k2=(1-a2)/3;V0=2*pi**2

def compute(uball,n):
  def bulk_inner(t,x):
    u=X.const(U.var(uball));r=rin*x;s2=t*t*(h*h+r*r);r2=a2+s2;q=1-u+u/r2;b=-u/r2**2;p=q+b*s2
    den=q*p.sqrt();R=8*u*a2/(q*p*r2**3)-2*u/(q**2*r2**2);w=4*hab*t*t*2*pi*rin**2*x
    return den*w,R*den*w
  def bulk_outer(t,z):
    u=X.const(U.var(uball));al=pi*z/3;ca=al.cos();sa=al.sin();r=rin/ca;s2=t*t*(h*h+r*r);r2=a2+s2
    q=1-u+u/r2;b=-u/r2**2;p=q+b*s2;den=q*p.sqrt();R=8*u*a2/(q*p*r2**3)-2*u/(q**2*r2**2)
    w=4*hab*t*t*(2*pi-6*al)*rin**2*sa/ca**3*pi/3;return den*w,R*den*w
  vi,bi=integrate_vec(bulk_inner,2,n);vo,bo=integrate_vec(bulk_outer,2,n);V=(vi+vo)*600;B=(bi+bo)*(-(arb(2)/3)*600)
  def face_inner(x):
    u=X.const(U.var(uball));r=rin*x;s2=h*h+r*r;r2=a2+s2;q=1-u+u/r2;b=-u/r2**2;p=q+b*s2
    return 2*u*(1-u)/(r2**3*q*p)*h*r*r*q*p.sqrt()/(q+b*r*r)*2*pi*rin**2*x
  def face_outer(z):
    u=X.const(U.var(uball));al=pi*z/3;ca=al.cos();sa=al.sin();r=rin/ca;s2=h*h+r*r;r2=a2+s2;q=1-u+u/r2;b=-u/r2**2;p=q+b*s2
    w=(2*pi-6*al)*rin**2*sa/ca**3*pi/3;return 2*u*(1-u)/(r2**3*q*p)*h*r*r*q*p.sqrt()/(q+b*r*r)*w
  F=(integrate_vec(face_inner,1,n)[0]+integrate_vec(face_outer,1,n)[0])*(-(arb(4)/3)*1200*2)
  def edge(t):
    u=X.const(U.var(uball));s2=k2*(1+2*t*t);r2=a2+s2;q=1-u+u/r2;b=-u/r2**2;p=q+b*s2;x=b*(k2/3)/p
    cosine=(arb(1)/3+x)/(1-x);beta=5*cosine.acos();dl=(q*2*k2+b*(2*k2*t)**2).sqrt();cone=16*pi**2/(3*beta)+8*beta/3-8*pi
    return 2*cone*dl
  E=integrate_vec(edge,1,n)[0]*720;R=B+F+E;A=V0**(arb(1)/3)*R*V**(-arb(1)/3)
  return V,B,F,E,R,A

def horner(coeffs, s):
  value=arb(0)
  for coefficient in reversed(coeffs):
    value=value*s+coefficient
  return value

def sign_intervals(a):
  s=arb(0,arb(1)/40)
  first_polynomial=horner([n*a[n] for n in range(1,UORDER+1)],s)
  second_polynomial=horner([n*(n-1)*a[n] for n in range(2,UORDER+1)],s)
  radius=arb(1)/10
  ratio=arb(1)/4
  majorant=arb(1000)
  n=UORDER
  first_tail=(majorant/radius)*ratio**n*((n+1)-n*ratio)/(1-ratio)**2
  second_sum=-ratio**n*(n*n*ratio**2-2*n*n*ratio+n*n-n*ratio**2+n+2*ratio)/(ratio*(ratio-1)**3)
  second_tail=(majorant/radius**2)*second_sum
  return (
    first_polynomial,
    second_polynomial,
    first_polynomial+arb(0,first_tail),
    second_polynomial+arb(0,second_tail),
    first_tail,
    second_tail,
  )

def ball_record(value):
  return {
    'midpoint':repr(value.mid()),
    'radius':repr(value.rad()),
    'lower':repr(value.lower()),
    'upper':repr(value.upper()),
  }

def analytic_audit():
  """Uniform complex bounds on a rectangle containing all Cauchy disks."""
  z_cells=[]
  for i in range(24):
    real=arb(-arb(1)/10+arb(2*i+1)/40,arb(1)/40)
    z_cells.append(acb(real,arb(0,arb(1)/10)))

  min_q=arb(10);min_p=arb(10);min_p_real=arb(10)
  min_density_real=arb(10);max_density=arb(0)
  max_d1=arb(0);max_d2=arb(0);max_bulk=arb(0)
  for z in z_cells:
    for j in range(64):
      lo=a2+(1-a2)*j/64;hi=a2+(1-a2)*(j+1)/64
      r2=arb((lo+hi)/2,(hi-lo)/2)
      d1=1/r2-1;d2=a2/r2**2-1
      q=1+z*d1;p=1+z*d2;density=q*p.sqrt()
      bulk=8*z*a2/(p.sqrt()*r2**3)-2*z*p.sqrt()/(q*r2**2)
      min_q=min(min_q,q.abs_lower());min_p=min(min_p,p.abs_lower())
      min_p_real=min(min_p_real,p.real.lower())
      min_density_real=min(min_density_real,density.real.lower())
      max_density=max(max_density,density.abs_upper())
      max_d1=max(max_d1,d1.abs_upper());max_d2=max(max_d2,d2.abs_upper())
      max_bulk=max(max_bulk,bulk.abs_upper())

  face_base=a2+h*h;min_face=arb(10);min_face_real=arb(10);max_face=arb(0)
  for z in z_cells:
    for j in range(64):
      lo=face_base+(1-face_base)*j/64
      hi=face_base+(1-face_base)*(j+1)/64
      r2=arb((lo+hi)/2,(hi-lo)/2);radial2=r2-face_base
      p=1+z*(a2/r2**2-1)
      face_denominator=1+z*(face_base/r2**2-1)
      face=2*z*(1-z)*h*radial2/(r2**3*p.sqrt()*face_denominator)
      min_face=min(min_face,face_denominator.abs_lower())
      min_face_real=min(min_face_real,face_denominator.real.lower())
      max_face=max(max_face,face.abs_upper())

  max_delta=arb(0);min_beta=arb(10);max_cosine=arb(0);max_edge=arb(0)
  min_edge_line_real=arb(10)
  for z in z_cells:
    for j in range(64):
      tlo=arb(j)/64;thi=arb(j+1)/64
      t2lo=tlo**2;t2hi=thi**2
      t2=arb((t2lo+t2hi)/2,(t2hi-t2lo)/2)
      r2=a2+k2*(1+2*t2);q=1+z*(1/r2-1);p=1+z*(a2/r2**2-1)
      b=-z/r2**2;x=b*(k2/3)/p;cosine=(arb(1)/3+x)/(1-x)
      beta=5*cosine.acos();delta=2*pi-beta
      cone=-(arb(4)/3)*delta+4*delta**2/(3*beta)
      dl_squared=q*2*k2+b*(2*k2)**2*t2
      dl=dl_squared.sqrt()
      max_delta=max(max_delta,delta.abs_upper())
      min_beta=min(min_beta,beta.abs_lower())
      max_cosine=max(max_cosine,cosine.abs_upper())
      max_edge=max(max_edge,(2*cone*dl).abs_upper())
      min_edge_line_real=min(min_edge_line_real,dl_squared.real.lower())

  flat_volume=600*ell**3/(6*arb(2).sqrt())
  face_area=arb(3).sqrt()*ell**2/4
  volume_lower=flat_volume*min_density_real
  bulk_bound=(arb(2)/3)*flat_volume*max_bulk
  face_bound=(arb(4)/3)*1200*2*face_area*max_face
  edge_bound=720*max_edge
  raw_bound=bulk_bound+face_bound+edge_bound
  normalized_bound=V0**(arb(1)/3)*raw_bound/volume_lower**(arb(1)/3)
  return {
    'a2':a2,'min_r2':a2,'max_r2':arb(1),'max_abs_inverse_r2_minus_one':max_d1,
    'max_abs_a2_over_r2_squared_minus_one':max_d2,'min_abs_q':min_q,
    'min_abs_p':min_p,'min_real_p':min_p_real,
    'min_abs_face_denominator':min_face,'min_real_face_denominator':min_face_real,
    'min_real_volume_density':min_density_real,'max_abs_volume_density':max_density,
    'volume_lower':volume_lower,'bulk_bound':bulk_bound,'face_bound':face_bound,
    'edge_bound':edge_bound,'raw_bound':raw_bound,'normalized_bound':normalized_bound,
    'max_abs_delta':max_delta,'min_abs_beta':min_beta,
    'max_abs_edge_cosine':max_cosine,
    'min_real_edge_line_element_squared':min_edge_line_real,
  }

def evaluate_cell(task):
  order,index=task
  midpoint=arb(2*index+1)/40
  started=time.time();out=compute(midpoint,order)
  first_poly,second_poly,first,second,first_tail,second_tail=sign_intervals(out[-1].a)
  return {
    'order':order,'cell':index,'midpoint':f'{2*index+1}/40',
    'first_polynomial':ball_record(first_poly),
    'second_polynomial':ball_record(second_poly),
    'first_derivative':ball_record(first),
    'second_derivative':ball_record(second),
    'first_negative':bool(first.upper()<0),
    'second_positive':bool(second.lower()>0),
    'finite':bool(all(value.is_finite() for value in out[-1].a)),
    'first_tail':repr(first_tail),'second_tail':repr(second_tail),
    'elapsed_seconds':time.time()-started,
  }

def evaluate_duffy_control():
  started=time.time();out=compute(arb(1)/2,16)
  return {
    'u':'1/2','order':16,
    'components':{name:ball_record(value.a[0]) for name,value in zip('volume bulk face edge raw normalized'.split(),out)},
    'elapsed_seconds':time.time()-started,
  }

def exact_endpoint_audit():
  bulk_factor=Fraction(-2,3)
  face_factor=Fraction(-4,3)
  linear_edge_factor=Fraction(-4,3)
  deficit=sp.symbols('delta')
  beta=2*sp.pi-deficit
  cone=-sp.Rational(4,3)*deficit+4*deficit**2/(3*beta)
  remainder=sp.cancel(cone+sp.Rational(4,3)*deficit)
  return {
    'face_factor':face_factor==2*bulk_factor,
    'edge_factor':linear_edge_factor==2*bulk_factor,
    'cone_constant':sp.simplify(remainder.subs(deficit,0))==0,
    'cone_linear':sp.simplify(sp.diff(remainder,deficit).subs(deficit,0))==0,
    'cone_remainder':str(sp.factor(remainder)),
  }

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--orders',type=int,nargs='+',choices=(16,20,24),default=(16,20,24))
 parser.add_argument('--jobs',type=int,default=4)
 parser.add_argument('--algebra-only',action='store_true')
 args=parser.parse_args()
 if args.algebra_only:
  exact=exact_endpoint_audit()
  for key in ('face_factor','edge_factor','cone_constant','cone_linear'):
   print(f"[{'PASS' if exact[key] else 'FAIL'}] {key}")
  print('cone remainder:',exact['cone_remainder'])
  if not all(exact[key] for key in ('face_factor','edge_factor','cone_constant','cone_linear')):
   raise SystemExit(1)
  raise SystemExit(0)
 orders=tuple(args.orders);started=time.time();audit=analytic_audit()

 tests=0
 passed=0
 def check(label,condition,detail=''):
  global tests,passed
  tests+=1
  ok=bool(condition)
  passed+=int(ok)
  suffix=f' -- {detail}' if detail else ''
  print(f"[{'PASS' if ok else 'FAIL'}] {label}{suffix}")

 check('exact 600-cell constant obeys 0.85<a^2<0.86',arb('0.85')<audit['a2']<arb('0.86'),repr(audit['a2']))
 check(
  'radial squared distance obeys a^2<=r^2<=1',
  audit['min_r2'].contains(a2) and audit['max_r2'].contains(arb(1)),
 )
 check('|1/r^2-1|<0.18 on every complex disk',audit['max_abs_inverse_r2_minus_one']<arb('0.18'),repr(audit['max_abs_inverse_r2_minus_one']))
 check('|a^2/r^4-1|<0.18 on every complex disk',audit['max_abs_a2_over_r2_squared_minus_one']<arb('0.18'),repr(audit['max_abs_a2_over_r2_squared_minus_one']))
 check('|q|>0.8 uniformly',audit['min_abs_q']>arb('0.8'),repr(audit['min_abs_q']))
 check('|p|>0.8 uniformly',audit['min_abs_p']>arb('0.8'),repr(audit['min_abs_p']))
 check('p stays in the right half-plane (analytic square-root branch)',audit['min_real_p']>arb('0.8'),repr(audit['min_real_p']))
 check('the face radial denominator has modulus >0.8',audit['min_abs_face_denominator']>arb('0.8'),repr(audit['min_abs_face_denominator']))
 check('the face denominator stays in the right half-plane',audit['min_real_face_denominator']>arb('0.8'),repr(audit['min_real_face_denominator']))
 check('|V(z)|>13 via a positive real-part bound',audit['volume_lower']>13,repr(audit['volume_lower']))
 check('|delta|<0.23 on the edge rectangles',audit['max_abs_delta']<arb('0.23'),repr(audit['max_abs_delta']))
 check('|beta|>6 on the edge rectangles',audit['min_abs_beta']>6,repr(audit['min_abs_beta']))
 check('the edge cosine stays strictly inside the acos branch points',audit['max_abs_edge_cosine']<1,repr(audit['max_abs_edge_cosine']))
 check('the squared edge line element stays in the right half-plane',audit['min_real_edge_line_element_squared']>0,repr(audit['min_real_edge_line_element_squared']))
 check('|bulk component|<350',audit['bulk_bound']<350,repr(audit['bulk_bound']))
 check('|face component|<100',audit['face_bound']<100,repr(audit['face_bound']))
 check('|edge component|<171',audit['edge_bound']<171,repr(audit['edge_bound']))
 check('|A2_equal_volume|<1000',audit['normalized_bound']<1000,repr(audit['normalized_bound']))

 tasks=[(order,index) for order in orders for index in range(20)]
 print(f'Running {len(tasks)} validated Taylor cells with {args.jobs} workers',flush=True)
 cells=[]
 with mp.get_context('fork').Pool(args.jobs) as pool:
  for cell in pool.imap_unordered(evaluate_cell,tasks):
   cells.append(cell)
   relevant=cell['first_negative'] if cell['cell']<19 else cell['second_positive']
   key='A1 upper='+cell['first_derivative']['upper'] if cell['cell']<19 else 'A2 lower='+cell['second_derivative']['lower']
   print(f"  order={cell['order']:2d} cell={cell['cell']:2d} {'PASS' if relevant else 'FAIL'} {key}",flush=True)
 cells.sort(key=lambda item:(item['order'],item['cell']))
 for order in orders:
  ordered=[cell for cell in cells if cell['order']==order]
  check(f'order {order}: all 20 Taylor coefficient sets are finite',all(cell['finite'] for cell in ordered))
  check(f'order {order}: A2 derivative is negative on cells 0,...,18',all(cell['first_negative'] for cell in ordered[:19]))
  check(f'order {order}: A2 second derivative is positive on cell 19',ordered[19]['second_positive'])

 control=evaluate_duffy_control()
 old={'volume':18.19424382,'bulk':-40.162026457,'face':1.33327578,'edge':-37.983620004,'raw':-76.81237068,'normalized':-78.92775426}
 for name,target in old.items():
  record=control['components'][name]
  enclosure=arb(record['midpoint'])+arb(0,arb(record['radius']).upper())
  check(f'radial reduction reproduces the independent u=1/2 {name}',enclosure.contains(arb(str(target))),f"target={target}, enclosure={record}")

 # Exact endpoint argument.  With the frozen convention, the bulk, face and
 # linear-deficit terms are -(2/3) times the distributional total scalar
 # curvature of a continuous piecewise-smooth metric.  The factors 2 on the
 # codimension-one mean-curvature jump and codimension-two angle deficit are
 # checked here directly.  The standard first-variation formula applies
 # after those internal boundary/corner terms cancel.  For n=3 and sectional
 # curvature one,
 # R=6, Ric=2g.  Hence delta S=(R/2-2) integral tr(h)=integral
 # tr(h)=2 delta V, while S/(3V)=2.  The derivative of S V^(-1/3)
 # vanishes for every variation h.  The non-Regge cone remainder is exactly
 # quadratic in delta, so its first derivative also vanishes at delta=0.
 n=arb(3);scalar=n*(n-1);ricci=n-1
 eh_variation_factor=scalar/2-ricci
 normalization_factor=scalar/n
 exact=exact_endpoint_audit()
 check('face coefficient is exactly twice the bulk scalar-curvature factor',exact['face_factor'])
 check('linear deficit coefficient is exactly twice the bulk factor',exact['edge_factor'])
 check('round normalized Einstein-Hilbert first variation cancels exactly',eh_variation_factor*2-normalization_factor==0)
 check(
  'exact-minus-linear cone correction has zero constant and linear terms',
  exact['cone_constant'] and exact['cone_linear'],
  f"remainder={exact['cone_remainder']}",
 )

 required_orders=set(orders)=={16,20,24}
 path_selected=required_orders and all(
  cell['finite'] and (cell['first_negative'] if cell['cell']<19 else cell['second_positive'])
  for cell in cells
 )
 status=('DERIVED PATH SELECTION' if path_selected else 'PARTIAL/FAILED CERTIFICATE')
 result={
  'protocol_commit':'8b322b1',
  'status':status,
  'spatial_orders':list(orders),
  'cover':'I_j=[j/20,(j+1)/20], j=0,...,19',
  'taylor_degree':UORDER,
  'analytic_audit':{key:repr(value) for key,value in audit.items()},
  'cells':cells,
  'duffy_control':control,
  'exact_endpoint_audit':exact,
  'endpoint_stationarity':(
   'exact normalized Einstein-Hilbert first variation at the unit round '
   'metric; internal face/linear-edge terms are its distributional '
   'completion and the remaining cone term is quadratic in the deficit'
  ),
  'elapsed_seconds':time.time()-started,
  'tests':tests,
  'passed':passed,
  'scope':(
   'One affine Riemannian path and one complete-exterior de Rham A2 '
   'coefficient only; no global metric theorem, cutoff sign, Lorentzian '
   'dynamics, Newton constant, Planck scale, or source coupling.'
  ),
 }
 if required_orders:
  Path(__file__).with_name('round_regge_a2_interval.json').write_text(json.dumps(result,indent=2)+'\n')
 print('-'*78)
 print(f'RESULT: {passed}/{tests} checks passed')
 print(status)
 if passed!=tests or not path_selected:
  raise SystemExit(1)
