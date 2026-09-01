#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

enum {
  RNDN = MPFR_RNDN,
  RNDZ = MPFR_RNDZ,
  RNDU = MPFR_RNDU,
  RNDD = MPFR_RNDD
};
enum { PREC = 256, MAX_VALUES = 400 };

typedef struct { mpfr_t lo, hi; } IV;

static void iv_init(IV *x) { mpfr_init2(x->lo, PREC); mpfr_init2(x->hi, PREC); }
static void iv_clear(IV *x) { mpfr_clear(x->lo); mpfr_clear(x->hi); }
static void iv_zero(IV *x) { mpfr_set_ui(x->lo, 0, RNDN); mpfr_set_ui(x->hi, 0, RNDN); }
static void iv_ui(IV *x, unsigned long n) { mpfr_set_ui(x->lo, n, RNDN); mpfr_set_ui(x->hi, n, RNDN); }
static void iv_str(IV *x, const char *s) {
  if (mpfr_set_str(x->lo, s, 10, RNDD) || mpfr_set_str(x->hi, s, 10, RNDU)) {
    fprintf(stderr, "bad MPFR decimal: %s\n", s); exit(2);
  }
}
static void iv_copy(IV *z, const IV *x) { mpfr_set(z->lo, x->lo, RNDD); mpfr_set(z->hi, x->hi, RNDU); }
static void iv_add(IV *z, const IV *x, const IV *y) { mpfr_add(z->lo,x->lo,y->lo,RNDD); mpfr_add(z->hi,x->hi,y->hi,RNDU); }
static void iv_sub(IV *z, const IV *x, const IV *y) { mpfr_sub(z->lo,x->lo,y->hi,RNDD); mpfr_sub(z->hi,x->hi,y->lo,RNDU); }
static void iv_mul(IV *z, const IV *x, const IV *y) {
  if (mpfr_cmp_ui(x->lo,0)<0 || mpfr_cmp_ui(y->lo,0)<0) { fprintf(stderr,"negative interval multiply\n"); exit(3); }
  mpfr_mul(z->lo,x->lo,y->lo,RNDD); mpfr_mul(z->hi,x->hi,y->hi,RNDU);
}
static void iv_div(IV *z, const IV *x, const IV *y) {
  if (mpfr_cmp_ui(x->lo,0)<0 || mpfr_cmp_ui(y->lo,0)<=0) { fprintf(stderr,"bad interval divide\n"); exit(3); }
  mpfr_div(z->lo,x->lo,y->hi,RNDD); mpfr_div(z->hi,x->hi,y->lo,RNDU);
}
static void iv_sqrt(IV *z, const IV *x) {
  if (mpfr_cmp_ui(x->lo,0)<0) { fprintf(stderr,"negative interval sqrt\n"); exit(3); }
  mpfr_sqrt(z->lo,x->lo,RNDD); mpfr_sqrt(z->hi,x->hi,RNDU);
}
static void iv_pow_gt1(IV *z, const IV *x, const IV *p) {
  mpfr_pow(z->lo,x->lo,p->lo,RNDD); mpfr_pow(z->hi,x->hi,p->hi,RNDU);
}
static void iv_scale_ui(IV *z, const IV *x, unsigned long n) {
  mpfr_t m; mpfr_init2(m,PREC); mpfr_set_ui(m,n,RNDN);
  mpfr_mul(z->lo,x->lo,m,RNDD); mpfr_mul(z->hi,x->hi,m,RNDU); mpfr_clear(m);
}
static double print_lower(mpfr_srcptr x) {
  double value=mpfr_get_d(x,RNDD);
  if(!isfinite(value)){fprintf(stderr,"non-finite lower endpoint\n");exit(3);}
  return nextafter(value,-INFINITY);
}
static double print_upper(mpfr_srcptr x) {
  double value=mpfr_get_d(x,RNDU);
  if(!isfinite(value)){fprintf(stderr,"non-finite upper endpoint\n");exit(3);}
  return nextafter(value,INFINITY);
}

typedef struct { IV x, logx; unsigned long count; int zero; } Value;
typedef struct { Value v[MAX_VALUES]; int m; unsigned long n; } Group;

static void read_groups(const char *path, Group g[2]) {
  FILE *f=fopen(path,"r"); if(!f){perror(path);exit(2);} char line[256];
  while(fgets(line,sizeof line,f)) {
    char *a=strtok(line,","), *b=strtok(NULL,","), *c=strtok(NULL,",\r\n");
    if(!a||!b||!c) { fprintf(stderr,"bad count line\n"); exit(2); }
    int k=atoi(a); if(k<0||k>1||g[k].m>=MAX_VALUES) exit(2);
    Value *v=&g[k].v[g[k].m++]; iv_init(&v->x); iv_init(&v->logx);
    iv_str(&v->x,b); v->count=strtoul(c,NULL,10);
    v->zero=mpfr_zero_p(v->x.lo) && mpfr_zero_p(v->x.hi);
    if(!v->zero && mpfr_cmp_ui(v->x.lo,1)<0) {
      fprintf(stderr,"positive spend below one: %s\n",b); exit(2);
    }
    if(v->zero) iv_zero(&v->logx); else {
      mpfr_log(v->logx.lo,v->x.lo,RNDD); mpfr_log(v->logx.hi,v->x.hi,RNDU);
    }
    g[k].n += v->count;
  }
  fclose(f);
}

static void group_mean_power(IV *out, const Group *g, const IV *p) {
  IV sum,pow,term,den; iv_init(&sum);iv_init(&pow);iv_init(&term);iv_init(&den);iv_zero(&sum);
  for(int i=0;i<g->m;i++) if(!g->v[i].zero) {
    iv_pow_gt1(&pow,&g->v[i].x,p); iv_scale_ui(&term,&pow,g->v[i].count); iv_add(&sum,&sum,&term);
  }
  iv_ui(&den,g->n); iv_div(out,&sum,&den);
  iv_clear(&sum);iv_clear(&pow);iv_clear(&term);iv_clear(&den);
}

static void group_envelope_means(IV *mh, IV *mh1, IV *mhh1, const Group *g, const IV *right) {
  IV sh,sh1,shh1,h,h1,hh1,term,den; IV *all[]={&sh,&sh1,&shh1,&h,&h1,&hh1,&term,&den};
  for(int i=0;i<8;i++) iv_init(all[i]); iv_zero(&sh);iv_zero(&sh1);iv_zero(&shh1);
  for(int i=0;i<g->m;i++) if(!g->v[i].zero) {
    iv_pow_gt1(&h,&g->v[i].x,right); iv_mul(&h1,&h,&g->v[i].logx); iv_mul(&hh1,&h,&h1);
    iv_scale_ui(&term,&h,g->v[i].count); iv_add(&sh,&sh,&term);
    iv_scale_ui(&term,&h1,g->v[i].count); iv_add(&sh1,&sh1,&term);
    iv_scale_ui(&term,&hh1,g->v[i].count); iv_add(&shh1,&shh1,&term);
  }
  iv_ui(&den,g->n); iv_div(mh,&sh,&den);iv_div(mh1,&sh1,&den);iv_div(mhh1,&shh1,&den);
  for(int i=0;i<8;i++) iv_clear(all[i]);
}

typedef struct { IV variance, raw2, lower, upper; } Node;
static void node_init(Node *x){iv_init(&x->variance);iv_init(&x->raw2);iv_init(&x->lower);iv_init(&x->upper);}
static void node_clear(Node *x){iv_clear(&x->variance);iv_clear(&x->raw2);iv_clear(&x->lower);iv_clear(&x->upper);}

static void make_rational(IV *out, unsigned long num, unsigned long den) {
  IV a,b;iv_init(&a);iv_init(&b);iv_ui(&a,num);iv_ui(&b,den);iv_div(out,&a,&b);iv_clear(&a);iv_clear(&b);
}

static void compute_node(Node *out, const Group g[2], const IV *p, const IV *wx, const IV *wy,
                         const IV *neff, const IV *q) {
  IV two,p2,mx,my,e2x,e2y,mx2,my2,vx,vy,t1,t2,v,d,ratio,root,hw;
  IV *all[]={&two,&p2,&mx,&my,&e2x,&e2y,&mx2,&my2,&vx,&vy,&t1,&t2,&v,&d,&ratio,&root,&hw};
  for(int i=0;i<17;i++)iv_init(all[i]);iv_ui(&two,2);iv_mul(&p2,p,&two);
  group_mean_power(&mx,&g[0],p);group_mean_power(&my,&g[1],p);
  group_mean_power(&e2x,&g[0],&p2);group_mean_power(&e2y,&g[1],&p2);
  iv_mul(&mx2,&mx,&mx);iv_mul(&my2,&my,&my);iv_sub(&vx,&e2x,&mx2);iv_sub(&vy,&e2y,&my2);
  iv_mul(&t1,wx,&vx);iv_mul(&t2,wy,&vy);iv_add(&v,&t1,&t2);iv_copy(&out->variance,&v);
  iv_mul(&t1,wx,&e2x);iv_mul(&t2,wy,&e2y);iv_add(&out->raw2,&t1,&t2);
  iv_sub(&d,&my,&mx);iv_div(&ratio,&v,neff);iv_sqrt(&root,&ratio);iv_mul(&hw,q,&root);
  iv_sub(&out->lower,&d,&hw);iv_add(&out->upper,&d,&hw);
  for(int i=0;i<17;i++)iv_clear(all[i]);
}

int main(int argc,char **argv){
  if(argc!=2){fprintf(stderr,"usage: %s counts.csv\n",argv[0]);return 2;}
  Group g[2]={{0}};read_groups(argv[1],g); if(g[0].n!=21306||g[1].n!=21307){fprintf(stderr,"bad n\n");return 2;}
  const int nc=4000; Node *nodes=calloc(nc+1,sizeof(Node));
  double *cell_lower=calloc(nc,sizeof(double)), *cell_upper=calloc(nc,sizeof(double));
  if(!nodes||!cell_lower||!cell_upper)return 2;
  IV wx,wy,neff,q,p,h,target,tmp1,tmp2,tmp3; IV *base[]={&wx,&wy,&neff,&q,&p,&h,&target,&tmp1,&tmp2,&tmp3};
  for(int i=0;i<10;i++)iv_init(base[i]);
  make_rational(&wx,g[1].n,g[0].n+g[1].n); make_rational(&wy,g[0].n,g[0].n+g[1].n);
  make_rational(&tmp1,g[0].n*g[1].n,g[0].n+g[1].n);iv_copy(&neff,&tmp1);
  iv_str(&q,"2.22823852982312");make_rational(&h,1,6400);iv_str(&target,"0.1");
  for(int k=0;k<=nc;k++) { node_init(&nodes[k]); make_rational(&p,800+k,3200); compute_node(&nodes[k],g,&p,&wx,&wy,&neff,&q); }

  IV min_ratio,max_excess,min_lower,min_vlo;
  iv_init(&min_ratio);iv_init(&max_excess);iv_init(&min_lower);iv_init(&min_vlo);
  iv_str(&min_ratio,"1e100");iv_zero(&max_excess);iv_str(&min_lower,"1e100");iv_str(&min_vlo,"1e100");
  int positive=0;
  for(int k=0;k<nc;k++) {
    IV right,hx,h1x,hh1x,hy,h1y,hh1y,px,py,lcx,lcy,lv,la,ld,vlo,aup,ratio,lw,lb,excess,celllo,cellup,den,sq;
    IV *z[]={&right,&hx,&h1x,&hh1x,&hy,&h1y,&hh1y,&px,&py,&lcx,&lcy,&lv,&la,&ld,&vlo,&aup,&ratio,&lw,&lb,&excess,&celllo,&cellup,&den,&sq};
    for(int i=0;i<24;i++)iv_init(z[i]); make_rational(&right,801+k,3200);
    group_envelope_means(&hx,&h1x,&hh1x,&g[0],&right);group_envelope_means(&hy,&h1y,&hh1y,&g[1],&right);
    iv_mul(&tmp1,&hx,&h1x);iv_add(&tmp2,&hh1x,&tmp1);iv_scale_ui(&lcx,&tmp2,2);
    iv_mul(&tmp1,&hy,&h1y);iv_add(&tmp2,&hh1y,&tmp1);iv_scale_ui(&lcy,&tmp2,2);
    iv_mul(&tmp1,&wx,&lcx);iv_mul(&tmp2,&wy,&lcy);iv_add(&lv,&tmp1,&tmp2);
    iv_mul(&tmp1,&wx,&hh1x);iv_mul(&tmp2,&wy,&hh1y);iv_add(&tmp3,&tmp1,&tmp2);iv_scale_ui(&la,&tmp3,2);
    iv_add(&ld,&h1x,&h1y);
    mpfr_set(vlo.lo, mpfr_cmp(nodes[k].variance.lo,nodes[k+1].variance.lo)<=0?nodes[k].variance.lo:nodes[k+1].variance.lo, RNDD);
    mpfr_set(vlo.hi, vlo.lo, RNDU);iv_mul(&tmp1,&lv,&h);mpfr_sub(vlo.lo,vlo.lo,tmp1.hi,RNDD);mpfr_set(vlo.hi,vlo.lo,RNDU);
    mpfr_set(aup.hi, mpfr_cmp(nodes[k].raw2.hi,nodes[k+1].raw2.hi)>=0?nodes[k].raw2.hi:nodes[k+1].raw2.hi, RNDU);
    mpfr_set(aup.lo,aup.hi,RNDD);iv_mul(&tmp1,&la,&h);mpfr_add(aup.hi,aup.hi,tmp1.hi,RNDU);mpfr_set(aup.lo,aup.hi,RNDD);
    if(mpfr_cmp_ui(vlo.lo,0)<=0){fprintf(stderr,"nonpositive vlo cell %d\n",k);return 4;}
    mpfr_div(ratio.lo,vlo.lo,aup.hi,RNDD);mpfr_set(ratio.hi,ratio.lo,RNDU);
    if(mpfr_cmp(ratio.lo,min_ratio.lo)<0)mpfr_set(min_ratio.lo,ratio.lo,RNDD);
    if(mpfr_cmp(vlo.lo,min_vlo.lo)<0)mpfr_set(min_vlo.lo,vlo.lo,RNDD);

    /* upper L_w = q L_v / {2 sqrt(N_eff v_lower)} */
    mpfr_mul(den.lo,neff.lo,vlo.lo,RNDD);mpfr_sqrt(sq.lo,den.lo,RNDD);mpfr_set_ui(tmp1.lo,2,RNDN);mpfr_mul(den.lo,sq.lo,tmp1.lo,RNDD);
    mpfr_mul(tmp2.hi,q.hi,lv.hi,RNDU);mpfr_div(lw.hi,tmp2.hi,den.lo,RNDU);mpfr_add(lb.hi,ld.hi,lw.hi,RNDU);
    mpfr_set_ui(tmp1.hi,3,RNDN);mpfr_mul(tmp2.hi,lb.hi,h.hi,RNDU);mpfr_mul(excess.hi,tmp2.hi,tmp1.hi,RNDU);
    if(mpfr_cmp(excess.hi,max_excess.hi)>0)mpfr_set(max_excess.hi,excess.hi,RNDU);
    mpfr_set(celllo.lo,mpfr_cmp(nodes[k].lower.lo,nodes[k+1].lower.lo)<=0?nodes[k].lower.lo:nodes[k+1].lower.lo,RNDD);
    mpfr_mul(tmp1.hi,lb.hi,h.hi,RNDU);mpfr_sub(celllo.lo,celllo.lo,tmp1.hi,RNDD);
    mpfr_set(cellup.hi,mpfr_cmp(nodes[k].upper.hi,nodes[k+1].upper.hi)>=0?nodes[k].upper.hi:nodes[k+1].upper.hi,RNDU);
    mpfr_add(cellup.hi,cellup.hi,tmp1.hi,RNDU);
    if(mpfr_cmp(celllo.lo,min_lower.lo)<0)mpfr_set(min_lower.lo,celllo.lo,RNDD);
    if(mpfr_cmp_ui(celllo.lo,0)>0)positive++;
    cell_lower[k]=mpfr_get_d(celllo.lo,RNDD);cell_upper[k]=mpfr_get_d(cellup.hi,RNDU);
    for(int i=0;i<24;i++)iv_clear(z[i]);
  }
  printf("n_x=%lu\nn_y=%lu\nprecision_bits=%d\ncells=%d\n",g[0].n,g[1].n,PREC,nc);
  printf("certified_guard_ratio_lower=%.17g\n",print_lower(min_ratio.lo));
  printf("certified_variance_lower=%.17g\n",print_lower(min_vlo.lo));
  printf("enclosure_excess_upper=%.17g\n",print_upper(max_excess.hi));
  printf("accuracy_target_upper=%.17g\n",print_upper(target.hi));
  printf("accuracy_target_met=%s\n",mpfr_cmp(max_excess.hi,target.lo)<=0?"TRUE":"FALSE");
  printf("minimum_whole_cell_lower_bound=%.17g\n",print_lower(min_lower.lo));
  printf("positive_cells=%d\nouter_cells=%d\n",positive,nc-positive);
  const int anchor_node[]={0,800,1600,2400,3200,4000};
  for(int a=0;a<6;a++) {
    int k=anchor_node[a]; double lo,up;
    if(k==0){lo=cell_lower[0];up=cell_upper[0];}
    else if(k==nc){lo=cell_lower[nc-1];up=cell_upper[nc-1];}
    else {lo=fmin(cell_lower[k-1],cell_lower[k]);up=fmax(cell_upper[k-1],cell_upper[k]);}
    printf("anchor_%d_order=%.2f\nanchor_%d_lower=%.17g\nanchor_%d_upper=%.17g\n",
           a+1,0.25+0.25*a,a+1,nextafter(lo,-INFINITY),
           a+1,nextafter(up,INFINITY));
  }

  for(int k=0;k<=nc;k++)node_clear(&nodes[k]);free(nodes);free(cell_lower);free(cell_upper);
  iv_clear(&min_ratio);iv_clear(&max_excess);iv_clear(&min_lower);iv_clear(&min_vlo);
  for(int i=0;i<10;i++)iv_clear(base[i]);
  for(int k=0;k<2;k++)for(int i=0;i<g[k].m;i++){iv_clear(&g[k].v[i].x);iv_clear(&g[k].v[i].logx);}
  return 0;
}
