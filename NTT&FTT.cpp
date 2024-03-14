//FFT
/*
lengli_QAQ
Hope there are no bugs!!!
*/
#include <bits/stdc++.h>
#define fastio ios::sync_with_stdio(0); cin.tie(0); cout.tie(0)
#define endl '\n'
//#define x first
//#define y second

using namespace std;
typedef pair<int,int> PII;
typedef long long LL;

const double pi=acos(-1);
const int N=3000010;

int R[N],mx,B;
int lsn;

struct com
{
	double x,y;
	com(){};
	com(double x,double y) : x(x),y(y){};
	friend com operator + (const com &a,const com &b){return com(a.x+b.x,a.y+b.y);}
	friend com operator - (const com &a,const com &b){return com(a.x-b.x,a.y-b.y);}
	friend com operator * (const com &a,const com &b){return com(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
};

void FFT(com *a,int n,int inv)
{
	if(lsn!=n and (lsn=n))
		for(int i=0;i<n;i++) R[i]=(R[i>>1]>>1) | ((i&1)<<(B-1));
	for(int i=0;i<n;i++) if(i<R[i]) swap(a[i],a[R[i]]);
	for(int i=1;i<n;i<<=1)
	{
		com mi(cos(pi/i),sin(pi/i)*inv);
		for(int j=0;j<n;j+=(i<<1))
		{
			com x(1,0);
			for(int k=0;k<i;k++,x=x*mi)
			{
				com t1=a[j+k],t2=x*a[j+k+i];
				a[j+k]=t1+t2;
				a[j+k+i]=t1-t2;
			}
		}
	}
	
}

int n,m;
com a[N],b[N],Z(0,0);

signed main()
{
    fastio;
    cin>>n>>m;
    for(B=0,mx=1;mx<=n+m;mx<<=1,B++);
    
 	for(int j=0;j<=mx;j++) a[j]=b[j]=Z;
    for(int i=0;i<=n;i++) cin>>a[i].x;
    for(int i=0;i<=m;i++) cin>>b[i].x;
    
 	
 	FFT(a,mx,1);
 	FFT(b,mx,1);   
    for(int i=0;i<mx;i++) a[i]=a[i]*b[i];
	FFT(a,mx,-1);
	
	for(int i=0;i<=n+m;i++) cout<<(int)(a[i].x/mx+0.5)<<" ";
    
    return 0;
}


//NTT

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

const double pi=acos(-1);
const int N=3000010,M=998244353;
const int G=3,Gx=332748118;

int R[N],mx,B;
int lsn;

int qmi(int a,int b){
	int res=1;
	while(b){
		if(b&1) res=res*a%M;
		a=a*a%M;
		b>>=1;
	}
	return res;
}

void NTT(int a[],int n,int inv){
	if(lsn!=n and (lsn=n))
		for(int i=0;i<n;i++) R[i]=(R[i>>1]>>1) | ((i&1)<<(B-1));
	for(int i=0;i<n;i++) if(i<R[i]) swap(a[i],a[R[i]]);
	for(int i=1;i<n;i<<=1){
		int mi=qmi(inv==1 ? G : Gx,(M-1)/(i<<1));
		for(int j=0;j<n;j+=(i<<1)){
			LL x=1;
			for(int k=0;k<i;k++,x=x*mi%M){
				int t1=a[j+k]%M,t2=x*a[j+k+i]%M;
				a[j+k]=(t1+t2)%M;
				a[j+k+i]=(t1-t2+M)%M;
			}
		}
	}
}

string l,r;
int n,m;
int a[N],b[N];

signed main()
{
    fastio;
	cin>>n>>m;
    for(B=0,mx=1;mx<=n+m;mx<<=1,B++);
    
    for(int i=0;i<=n;i++) cin>>a[i],a[i]%=M;
    for(int i=0;i<=m;i++) cin>>b[i],b[i]%=M;
    
 	NTT(a,mx,1);
 	NTT(b,mx,1);   
    for(int i=0;i<mx;i++) a[i]=a[i]*b[i]%M;
	NTT(a,mx,-1);
	int inv=qmi(mx,M-2)%M;

	for(int i=0;i<=n+m;i++) cout<<(a[i]*inv)%M<<" ";
    
    return 0;
}

//FWT

namespace FWT{//idx:0->(1<<n)
	void FWT_OR(Mint *a,int op){
		for(int i=1;i<(1<<n);i*=2){//op=1,-1
			for(int p=i*2,j=0;j<(1<<n);j+=p){
				for(int k=0;k<i;k++){
					a[i+j+k]+=a[j+k]*op;
				}
			}
		}
	}
	void FWT_AND(Mint *a,int op){//op=1,-1
		for(int i=1;i<(1<<n);i*=2){
			for(int p=i*2,j=0;j<(1<<n);j+=p){
				for(int k=0;k<i;k++){
					a[j+k]+=a[i+j+k]*op;
				}
			}
		}
	}
	void FWT_XOR(Mint *a,int op){//op=1,inv(2)
		for(int i=1;i<(1<<n);i*=2){
			for(int p=i*2,j=0;j<(1<<n);j+=p){
				for(int k=0;k<i;k++){
					Mint x=a[j+k],y=a[i+j+k];
					a[j+k]=(x+y)*op;
					a[i+j+k]=(x-y)*op;
					
				}
			}
		}
	}
};
