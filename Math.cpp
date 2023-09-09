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



//线性筛质数与欧拉筛
#include <bits/stdc++.h>
//#define int long long
#define x first
#define y second
using namespace std;

typedef pair<int,int> PII;
typedef long long LL;
const int N=1010;

int n;
int prime[N],st[N],cnt;
int phi[N];

void ola()
{
    phi[1]=1;
    for(int i=2;i<=N;i++)
    {
        if(!st[i])
        {
            prime[cnt++]=i;
            st[i]=true;
            phi[i]=i-1;
        }
        for(int j=0;prime[j]*i<=N;j++)
        {
            st[prime[j]*i]=true;
            if(i%prime[j]==0)
            {
                phi[i*prime[j]]=prime[j]*phi[i];
                break;
            }
            phi[i*prime[j]]=(prime[j]-1)*phi[i];
        }
    }
}

signed main()
{
    //ios::sync_with_stdio(false);
    int T;
    cin>>T;
    int t=1;
    ola();
    while(T--)
    {
        cin>>n;
        LL res=0;
        for(int i=1;i<=n;i++) res+=phi[i]*2;
        cout<<t<<" "<<n<<" "<<res+1<<endl;
        t++;
    }
}

//快速幂

int qmi(int a,int b)
{
    int res=1;
    while(b)
    {
        if(b&1) res=(res*a)%M;
        a=(a*a)%M;
        b>>=1;
    }
    return res;
}

//矩阵快速幂

int A[4][4]={
    {0,1,0,0},
    {1,1,1,0},
    {0,0,1,1},
    {0,0,0,1}
};
int B[4][4]={
    {1,0,0,0},
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1}
};
int ini[4]={1,1,1,0};
int res[4]={0,0,0,0};

void mul(int (&a)[4][4],int (&b)[4][4])
{
    int C[4][4];
    memset(C,0,sizeof C);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            for(int k=0;k<4;k++)
            {
                C[i][j]=(C[i][j]+a[i][k]*b[k][j])%m;
            }
        }
    }
    memcpy(a,C,sizeof C);
}

void qmi(int b)
{
    while(b)
    {
        if(b&1) mul(B,A);
        mul(A,A);
        b>>=1;
    }
}

signed main()
{
    fastio;
    cin>>n>>m;
    if(n<=2)
    {
        if(n==1) cout<<1%m;
        else cout<<3%m;
    }
    else
    {
        qmi(n-1);
        
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                 res[i]=(res[i]+ini[j]*B[j][i])%m;
            }
        }
 
        
        cout<<(n*res[2]-res[3]+m)%m;
    }
}

//线性基

#include <cstring>
#include <iostream>
#include <algorithm>
#define int long long

using namespace std;
const int N=100010;

int n;
int a[N];
int p[N];

void add(int x)//获取线性基
{
    for(int i=62;i>=0;i--)
    {
        if(x>>i&1)
        {
            if(p[i])x^=p[i];
            else
            {
                p[i]=x;
                break;
            }
        }
    }
}

signed main()
{
    cin>>n;
    for(int i=0;i<n;i++) cin>>a[i];

    int res=0;
    for(int i=0;i<n;i++)
    {
        add(a[i]);
    }
    for(int i=62;i>=0;i--) 
    {
        if(p[i] and !(res>>i&1)) res^=p[i];
    }
    cout<<res;
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

//扩展欧几里得(ax+by==gcd(a,b))
namespace Exgcd{
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
