//字符串哈希
struct string_hash{
	int b=13331, P=1e9+7, p[N], rh[N],h[N],len;
	int inline get(int l, int r){
	    return (h[r] - (LL)h[l - 1] * p[r - l + 1] % P + P) % P;
	}
	void inline build(int n,string s) {
		p[0] = 1,len=n;
	    for(int i = 1; i <= n; i++){
	        p[i] = (LL)p[i - 1] * b % P;
	        h[i] = ((LL)h[i - 1] * b + s[i-1]) % P;
	    }
	    for(int i=1;i<=n;i++) rh[i]=((LL)rh[i - 1] * b + s[n-i]) % P;
	}
	bool isprs(int l,int r){
		return (h[r]-(LL)h[l-1]*p[r-l+1]%P+P)%P==(rh[len-l+1]-(LL)rh[len-r]*p[r-l+1]%P+P)%P;
	}
}a;

//KMP

vector<int> KMP(string a,string b){
    int n=a.size(),m=b.size();
	a="#"+a,b="#"+b;
	vector<int> q;
    vector<int> ne(m+1,0);
	for(int i=2,j=0;i<=m;i++){
		while(j and b[j+1]!=b[i]) j=ne[j];
		if(b[j+1]==b[i]) j++;
		ne[i]=j;
	}
	for(int i=1,j=0;i<=n;i++){
		while(j and a[i]!=b[j+1]) j=ne[j];
		if(a[i]==b[j+1]) j++;
		q[i]=j;
		if(j==m) j=ne[j];
	}
	return q;
}

//Z函数（从i开始的后缀与原串最长公共前缀）

z[1]=s.size()-1;
for(int i=2,l,r=0;i<s.size();i++){
	if(i<=r) z[i]=min(z[i-l+1],r-i+1);
	while(z[i]+i<s.size() and s[z[i]+1]==s[z[i]+i]) z[i]++;
	if(i+z[i]-1>r) r=i+z[i]-1,l=i;
}

//后缀数组（SA）
#include <bits/stdc++.h>
#define fastio ios::sync_with_stdio(0); cin.tie(0); cout.tie(0)
#define endl '\n'
//#define x first
//#define y second

using namespace std;
typedef pair<int,int> PII;
typedef long long LL;

const int N=1000010;

int n,m;
string s;
int sa[N], x[N], y[N], c[N], rk[N], height[N];

void get_sa()
{
    for (int i = 1; i <= n; i ++ ) c[x[i] = s[i]] ++ ;
    for (int i = 2; i <= m; i ++ ) c[i] += c[i - 1];
    for (int i = n; i; i -- ) sa[c[x[i]] -- ] = i;
    for (int k = 1; k <= n; k <<= 1)
    {
        int num = 0;
        for (int i = n - k + 1; i <= n; i ++ ) y[ ++ num] = i;
        for (int i = 1; i <= n; i ++ )
            if (sa[i] > k)
                y[ ++ num] = sa[i] - k;
        for (int i = 1; i <= m; i ++ ) c[i] = 0;
        for (int i = 1; i <= n; i ++ ) c[x[i]] ++ ;
        for (int i = 2; i <= m; i ++ ) c[i] += c[i - 1];
        for (int i = n; i; i -- ) sa[c[x[y[i]]] -- ] = y[i], y[i] = 0;
        swap(x, y);
        x[sa[1]] = 1, num = 1;
        for (int i = 2; i <= n; i ++ )
            x[sa[i]] = (y[sa[i]] == y[sa[i - 1]] && y[sa[i] + k] == y[sa[i - 1] + k]) ? num : ++ num;
        if (num == n) break;
        m = num;
    }
}

void get_height()
{
    for (int i = 1; i <= n; i ++ ) rk[sa[i]] = i;
    for (int i = 1, k = 0; i <= n; i ++ )
    {
        if (rk[i] == 1) continue;
        if (k) k -- ;
        int j = sa[rk[i] - 1];
        while (i + k <= n && j + k <= n && s[i + k] == s[j + k]) k ++ ;
        height[rk[i]] = k;
    }
}


void solve()
{
	cin>>s;
	s="#"+s;
	n=s.size()-1;
	m=125;
	get_sa();
	get_height();
	for(int i=1;i<=n;i++) cout<<sa[i]-1<<" ";
	cout<<"\n";
	for(int i=1;i<=n;i++) cout<<height[i]<<" ";
}

signed main()
{
    fastio;
    
    int T;
    T=1;
    while(T--) solve();
    
    return 0;
}

//Tire树
struct Trie {
  int nex[N * 32][2], idx;
  int cnt[N * 32];
  void insert(int x, int k) {
    int p = 0;
    for (int i = 31; i >= 0; i--) {
      int c = (x >> i) & 1;
      if (!nex[p][c]) nex[p][c] = ++idx;
      p = nex[p][c];
      cnt[p] += k;
    }
  }
 
  int find(int x) {
    int p = 0;
    int res = 0;
    for (int i = 31; i >= 0; i--) {
      int c = (x >> i) & 1;
      if (nex[p][c ^ 1] and cnt[nex[p][c ^ 1]]) {
        res += (1 << i);
        p = nex[p][c ^ 1];
      } else p = nex[p][c];
    }
    return res;
  }
} tr;