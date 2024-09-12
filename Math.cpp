//逆元预处理
void init_inv()
{
	inv[1] = 1;	
	for(int i = 2;i < N;i++) inv[i] = (M -M / i) * inv[M % i] %M;
}

//组合数（C）预处理

int aa[N],bb[N];
int inv[N];
 
void init(){
	inv[1] = 1;
	for(int i = 2;i < N;i++)  inv[i] = (M -  M / i) * inv[M % i] % M;
	aa[0]=1,bb[0]=1;
	for(int i=1;i<N;i++) aa[i]=(i*aa[i-1])%M;
	for(int i=1;i<N;i++) bb[i]=(inv[i]*bb[i-1])%M;
}
 
int C(int x, int y) {
	return x < y || y < 0 ? 0 : aa[x] * bb[y] % M * bb[x - y] % M;
}


//中国剩余定理excrt
namespace excrt{
	constexpr long long safe_mod(long long x, long long m) {
    x %= m;
    if (x < 0) x += m;
    return x;
	}
	std::pair<long long, long long> inv_gcd(long long a, long long b) {
	    a = safe_mod(a, b);
	    if (a == 0) return {b, 0};
	    long long s = b, t = a;
	    long long m0 = 0, m1 = 1;
	
	    while (t) {
	        long long u = s / t;
	        s -= t * u;
	        m0 -= m1 * u;  
	        auto tmp = s;
	        s = t;
	        t = tmp;
	        tmp = m0;
	        m0 = m1;
	        m1 = tmp;
	    }
	    if (m0 < 0) m0 += b / s;
	    return {s, m0};
	}
	pair<long long, long long> crt(const vector<long long>& r,const vector<long long>& m) {//a->val,m->mod
	    assert(r.size() == m.size());
	    int n = int(r.size());
	    long long r0 = 0, m0 = 1;
	    for (int i = 0; i < n; i++) {
	        assert(1 <= m[i]);
	        long long r1 = safe_mod(r[i], m[i]), m1 = m[i];
	        if (m0 < m1) {
	            std::swap(r0, r1);
	            std::swap(m0, m1);
	        }
	        if (m0 % m1 == 0) {
	            if (r0 % m1 != r1) return {0, 0};
	            continue;
	        }
	        long long g, im;
	        tie(g, im) = inv_gcd(m0, m1);
	        long long u1 = (m1 / g);
	        if ((r1 - r0) % g) return {0, 0};
	        long long x = (r1 - r0) / g % u1 * im % u1;
	        r0 += x * m0;
	        m0 *= u1;  // -> lcm(m0, m1)
	        if (r0 < 0) r0 += m0;
	    }
	    return {r0, m0};
	}
};

//线性筛质数

vector<int> init_prime(int n){
	vector<int> prime,st(n+1,0);
	for(int i=2;i<=n;i++){
		if(!st[i]) prime.push_back(i);
		for(int j=0;prime[j]*i<=n;j++){
			st[i*prime[j]]=1;
			if(i%prime[j]==0) break;
		}
	}
	return prime;
}

//线性处理欧拉函数

vector<int> init_phi(int n){
	vector<int> prime,st(n+1,0);
	vector<int> phi(n+1,0);
	phi[1]=1;
	for(int i=2;i<=n;i++){
		if(!st[i]) prime.push_back(i),phi[i]=i-1;
		for(int j=0;prime[j]*i<=n;j++){
			st[i*prime[j]]=1;
			if(i%prime[j]==0) {
				phi[i*prime[j]]=phi[i]*prime[j];
				break;
			}
			phi[i*prime[j]]=phi[i]*(prime[j]-1);
		}
	}
	return phi;
}

//欧拉降幂(欧拉定理)

struct ola_pow{
	int P,phi;
	void init(int p){
		P=p,phi=Phi(p);
	}
	int Phi(int n){
	    int res=n,m=n;
	    for(int i=2;i<=sqrt(m);i++){
	        if(m%i==0) res=res/i*(i-1);
	        while(m%i==0) m/=i;
	    }
	    if(m>1) res=res/m*(m-1);
	    return res;
	}
	int qmi(int a,int b,int p){//(a^b)%p
		int res=1;
		while(b){
			if(b&1) res=res*a%p;
			a=a*a%p;
			b>>=1;
		}
		return res%p;
	}
	int ola(int a,int b,int c){//(a^(b^c))%P
		bool flag=0;
		int k=1,d=(c>=10 ? 10 : c);
		for(int i=1;i<=d;i++) {
			k*=b;
			if(k>=phi) k%=phi,flag=1;
		}
		k*=qmi(b,c-d,phi);k%=phi;
		if(flag) k+=phi;
		return qmi(a%P,k,P);
	}
	int ola(int a,string s){//(a^s)%P
		bool flag=0;
		int res=0;
		for(int i=0;i<s.size();i++){
			int t=s[i]-'0';
			res*=10,res+=t;
			if(res>=phi) res%=phi,flag=1;
		}
		if(flag) res+=phi;
		return qmi(a,res,P);
	}
};

//矩阵快速幂

int h=2;
struct Matrix{
    vector<vector<int>> a = vector<vector<int>>(h,vector<int>(h));
    Matrix operator * (Matrix b){
        Matrix res;
        for(int i=0;i<h;i++)
            for(int j=0;j<h;j++)
                for(int k=0;k<h;k++) res.a[i][k]=(res.a[i][k]+(long long)a[i][j]*b.a[j][k])%M;
        return res;
    }
};
Matrix pow(Matrix a,long long b){
    Matrix res;
    for(int i=0;i<h;i++) res.a[i][i]=1;
    while(b){
        if(b&1) res=res*a;
        a=a*a;
        b>>=1; 
    }
    return res;
}

//线性基

struct Linear_basis{
	#define int long long
	std::vector<int> p;
	std::vector<int> b;
	int cnt=0;
	bool flag=0;//have zero ? 1 : 0
	void insert(int x){
		for(int i=62;i>=0;i--){
			if(x>>i&1){
				if(!p[i]){
					p[i]=x;
					return;
				}else x^=p[i];
			}
		}
		flag=1;
	}
	void init(std::vector<int> &a){
		p.clear();p.resize(63,0);
		flag=0;cnt=0;b.clear();
		for(auto x:a) insert(x);
		for(int i=0;i<=62;i++){
			for(int j=i-1;j>=0;j--)
				if(p[i]>>j&1) p[i]^=p[j];
			if(p[i]) cnt++,b.push_back(p[i]);
		}
	}
	int query_mi(){
		if(!flag) return 0;
		for(int i=0;i<=62;i++) 
			if(p[i]) return p[i];
	}
	int query_ma(){
		int res=0;
		for(int i=62;i>=0;i--)
			res=std::max(res,res^p[i]);
		return res;
	}
	int query(int k){//k>=1
		int res=0;
		k-=flag;
		if(!k) return 0;
		if(k>=(1ll<<cnt)) return -1;
		for(int i=0;i<cnt;i++)
			if(k>>i&1) res^=b[i];
		return res;
	}
	int query_rk(int x){
		int res=0;
		for(int i=cnt-1;i>=0;i--)
			if(x>=b[i]) res+=(1ll<<i),x^=b[i];
		return res+flag;
	}
}lb;

//质数检验

namespace miller{
	using u32=uint32_t;
	using u64=uint64_t;
	using u128=__uint128_t;
	constexpr u64 get_nr(u64 M){
	    u64 IV=2-M;
	    for(int i=0;i<5;++i){
	        IV*=2-M*IV;
	    }
	    return IV;
	}
	constexpr u64 mul(u64 x,u64 y,u64 IV,u64 M){
	    auto t=u128(x)*y;
	    u64 lo=t,hi=t>>64,res=(hi+M)-u64((u128(lo*IV)*M)>>64);
	    return res;
	}
	constexpr bool isprime(u64 x){
	    if(x<64){return (u64(1)<<x)&0x28208a20a08a28ac;}
	    if(x%2==0){return false;}
	    u64 d=x-1,IV=get_nr(x),R=(-x)%x,R2=(-u128(x))%x,nR=x-R;
	    int k=__builtin_ctzll(d);
	    d>>=k;
	    auto mr3=[&](u64 t1,u64 t2,u64 t3){
	        u64 r1=R,r2=R,r3=R;
	        t1=mul(t1,R2,IV,x),t2=mul(t2,R2,IV,x),t3=mul(t3,R2,IV,x);
	        for(u64 b=d;b;b>>=1){
	            if(b&1){
	                r1=mul(r1,t1,IV,x),r2=mul(r2,t2,IV,x),r3=mul(r3,t3,IV,x);
	            }
	            t1=mul(t1,t1,IV,x),t2=mul(t2,t2,IV,x),t3=mul(t3,t3,IV,x);
	        }
	        r1=std::min(r1,r1-x),r2=std::min(r2,r2-x),r3=std::min(r3,r3-x);
	        int res1=(r1==R)|(r1==nR),res2=(r2==R)|(r2==nR),res3=(r3==R)|(r3==nR);
	        for(int j=0;j<k-1;++j){
	            r1=mul(r1,r1,IV,x),r2=mul(r2,r2,IV,x),r3=mul(r3,r3,IV,x);
	            res1|=(std::min(r1,r1-x)==nR),res2|=(std::min(r2,r2-x)==nR),res3|=(std::min(r3,r3-x)==nR);
	        }
	        return res1&res2&res3;
	    };
	    auto mr4=[&](u64 t1,u64 t2,u64 t3,u64 t4){
	        u64 r1=R,r2=R,r3=R,r4=R;
	        t1=mul(t1,R2,IV,x),t2=mul(t2,R2,IV,x),t3=mul(t3,R2,IV,x),t4=mul(t4,R2,IV,x);
	        for(u64 b=d;b;b>>=1){
	            if(b&1){
	                r1=mul(r1,t1,IV,x),r2=mul(r2,t2,IV,x),r3=mul(r3,t3,IV,x),r4=mul(r4,t4,IV,x);
	            }
	            t1=mul(t1,t1,IV,x),t2=mul(t2,t2,IV,x),t3=mul(t3,t3,IV,x),t4=mul(t4,t4,IV,x);
	        }
	        r1=std::min(r1,r1-x),r2=std::min(r2,r2-x),r3=std::min(r3,r3-x),r4=std::min(r4,r4-x);
	        int res1=(r1==R)|(r1==nR),res2=(r2==R)|(r2==nR),res3=(r3==R)|(r3==nR),res4=(r4==R)|(r4==nR);
	        for(int j=0;j<k-1;++j){
	            r1=mul(r1,r1,IV,x),r2=mul(r2,r2,IV,x),r3=mul(r3,r3,IV,x),r4=mul(r4,r4,IV,x);
	            res1|=(std::min(r1,r1-x)==nR),res2|=(std::min(r2,r2-x)==nR),res3|=(std::min(r3,r3-x)==nR),res4|=(std::min(r4,r4-x)==nR);
	        }
	        return res1&res2&res3&res4;
	    };
	    if(x<(u64(1)<<32)){
	        return mr3(2,7,61);
	    }
	    return mr3(2,325,9375)&&mr4(28178,450775,9780504,1795265022);
	}
}


//超大质因数分解
/*
lengli_QAQ
Hope there are no bugs!!!
*/
#include <bits/stdc++.h>
#define fastio ios::sync_with_stdio(0); cin.tie(0); cout.tie(0)
#define endl '\n'
#define int long long
//#define x first
//#define y second

using namespace std;
typedef pair<int,int> PII;
typedef long long LL;

LL mul( LL a , LL b , LL md ) {
    return ( a * b - (LL)( (long double)a / md * b + 0.5 ) * md + md ) % md;
}
LL Pow( LL x , LL a , LL md ) {
    LL cur = x % md , ans = 1;
    while( a ) {
        if( a & 1 ) ans = mul( ans , cur , md );
        cur = mul( cur , cur , md ) , a >>= 1;
    }
    return ans;
}

const int ck[] = {2,3,5,7,11,13,17,19,23,29,31,37,41} , _l = 10;
bool miLLer( LL n ) {
    if( n == 1 ) return false;
    LL t = n - 1; int cn = 0;
    while( !( t & 1 ) ) t >>= 1 , ++ cn;
    for( int i = 0 ; i < _l ; ++ i ) {
        if( n == ck[i] ) return true;
        LL a = Pow( ck[i] , t , n ) , nex = a;
        for( int j = 1 ; j <= cn ; ++ j ) {
            nex = mul( a , a , n );
            if( nex == 1 && a != 1 && a != n - 1 ) return false;
            a = nex;
        }
        if( a != 1 ) return false;
    }
    return true;
}

inline LL f( LL x , LL c , LL md ) {
    return ( mul( x , x , md ) + c ) % md;
}

inline LL _rand(  ) {
    return (LL) rand() << 32 | rand();
}
inline LL _randw() {
    return (LL)rand() << 48 | (LL)rand() << 32 | rand() << 16 | rand();
}
inline LL _abs( LL x ) {
    return x > 0 ? x : -x;
}
inline LL gcd( LL a , LL b ) {
    return !b ? a : gcd( b , a % b );
}

inline LL poLLard_rho( LL n ) {
    LL s = 0 , t = 0 , c = _rand() % ( n - 1 ) + 1 , val = 1;
    for( int cir = 1 ; ; cir <<= 1 , s = t , val = 1 ) {
        for( int i = 0 ; i < cir ; ++ i ) {
            t = f( t , c , n ) , val = mul( val , _abs( t - s ) , n );
            if( i % 127 == 0 ) {
                LL g = gcd( val , n );
                if( g != 1 ) return g;
            }
        }
        LL g = gcd( val , n );
        if( g != 1 ) return g;
    }
}

vector<LL> divs;
inline void analyze( LL n ) {
    if( n == 1 ) return;
    if( miLLer( n ) ) { divs.push_back( n ); return; }
    LL d = n;
    while( d == n ) d = poLLard_rho( n );
    n /= d;
    analyze( n ) , analyze( d );
}

signed main()
{
    fastio;
    srand(time(0));


    
    return 0;
}

//高斯消元求线性方程组

int gauss(vector<vector<double>> &q){
	int n=q.size();
	assert(n==q[0].size()-1);
    int c,r;
    for(c=0,r=0;c<n;c++){
        int t=r;
        for(int i=r;i<n;i++)
           if(fabs(q[i][c])>fabs(q[t][c])) t=i;
        if(fabs(q[t][c])<eps) continue;
        for(int i=c;i<=n;i++) swap(q[t][i],q[r][i]);
        for(int i=n;i>=c;i--) q[r][i]/=q[r][c];
        for(int i=r+1;i<n;i++)
            if(fabs(q[i][c])>eps)
               for(int j=n;j>=c;j--)
                  q[i][j]-=q[r][j]*q[i][c];      
        r++;
    }
    if(r<n){
        for(int i=r;i<n;i++)
           if(fabs(q[i][n])>eps)
              return 2;
        return 1;
    }
    for(int i=n-1;i>=0;i--)
       for(int j=i+1;j<n;j++)
          q[i][n]-=q[i][j]*q[j][n];
    
    return 0;
}


//扩展欧几里得(ax+by==gcd(a,b))
namespace Exgcd{//通解x=x0+b/d,y=y0-a/d
	i128 x,y,a,b,d;
	#define chmax(a,b) (a>b ? a:b)
	#define chmin(a,b) (a>b ? b:a)
	void init(i128 aa,i128 bb){
		a=aa,b=bb;
		x=0,y=0,d=0;
	}
	i128 exgcd(i128 a,i128 b){
		if(!b){
			x=1,y=0;
			return a;
		}
		d=exgcd(b,a%b);
		i128 z=x;
		x=y;
		y=z-y*(a/b);
		return d;
	}
	void minz(){
		i128 t=0;
		if(x<0) t=(-x)/(b/d)+((-x)%(b/d)!=0);
		if(y<0) t=chmax(t,(-y)/(a/d)+((-y)%(a/d)!=0));
		x+=t*b/d;
		y+=t*a/d;
	}
	bool get(i128 c){
		x*=(c/d),y*=(c/d);
	    i128 dx=b/d,dy=a/d;
	    i128 t = 0;
	  	if (x < 0) t =(-x) / (b / d) + ((-x) % (b / d) != 0);
	  	if (y > 0) t =chmin(t, (y) / (a / d) + ((y) % (a / d) != 0));
	    x+=t*dx,y-=t*dy;
	    t=x/dx;
		x-=t*dx,y+=t*dy;
		if(x*a+y*b!=c){
			return 0;
		}
		return 1;
	}
}
