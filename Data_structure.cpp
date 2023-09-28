//ST表

struct ST{
	int a[N];
	int f[N][M];
	void init(){
		for(int len=0;len<M;len++){
			for(int i=1;i+(1<<len)-1<=n;i++){
				if(!len) f[i][len]=a[i];
				else f[i][len]=max(f[i][len-1],f[i+(1<<len-1)][len-1]);
			}
		}
	}
	int query(int l,int r){
		int k=log(r-l+1)/log(2);
		return max(f[l][k],f[r-(1<<k)+1][k]);
	}
}st;

//并查集

struct DSU{
    std::vector<int> p, sz,add;
    DSU(int n): p(n), sz(n, 1),add(n,0){
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x){//如果使用add，那么不能路径压缩！
        return x == p[x] ? x : p[x]=find(p[x]);
    }
    int sum(int x){
    	return x == p[x] ? add[x] : add[x]+sum(p[x]);
    }
    bool same(int x, int y){
        return find(x) == find(y);
    }
    bool merge(int x, int y){
        x = find(x);
        y = find(y);
        if (x == y) return 0;
        if(size(x) < size(y)) swap(x,y);
        sz[x] += sz[y];
        p[y] = x;
        add[y]-=add[x];
        return 1;
    }
    int size(int x){
        return sz[find(x)];
    }
};

//树状数组
struct BIT{
	int n=N-1;
	int a[N];
	int lowbit(int x){return x&-x;}
	void add(int x,int c) {for(int i=x;i<=n;i+=lowbit(i)) a[i]+=c;}
	
	void init(){
		
	};
	
	int sum(int x){
		LL res=0;
		for(int i=x;i;i-=lowbit(i)) res+=a[i];
		return res;
	}
	
}tr;

//珂朵莉树
struct Node
{
	int l,r;
	mutable int v;
	Node(int L,int R,int V) {l=L,r=R,v=V;}
	bool operator < (const Node &W) const {return l<W.l;}
};

set<Node> s;

auto split(int x)//分割区间
{
	auto it=s.lower_bound(Node(x,0,0));
	if(it!=s.end() and it->l == x) return it;
	it--;
	int l=it->l,r=it->r,v=it->v;
	s.erase(it);
	s.insert(Node(l,x-1,v));
	return s.insert(Node(x,r,v)).first;
}

void assign(int l,int r,int v)//区间平推赋值
{
	auto it2=split(r+1),it1=split(l);
	s.erase(it1,it2);
	s.insert(Node(l,r,v));
}

//线段树

#include <cstring>
#include <iostream>
#include <algorithm>

using namespace std;

typedef long long LL;
const int N=100010;

int n,m;
int a[N];

struct Node{
    int l,r;
    LL sum,add;
}tr[N*4];

void pushup(Node &u,Node &l,Node &r){
    u.sum=l.sum+r.sum;
}

void pushup(int u){
    pushup(tr[u],tr[u<<1],tr[u<<1|1]);
}

void pushdown(int u){
    if(tr[u].add){
        tr[u<<1].add+=tr[u].add,tr[u<<1].sum+=(LL)(tr[u<<1].r-tr[u<<1].l+1)*tr[u].add;
        tr[u<<1|1].add+=tr[u].add,tr[u<<1|1].sum+=(LL)(tr[u<<1|1].r-tr[u<<1|1].l+1)*tr[u].add;
        tr[u].add=0;
    }
}

void build(int u,int l,int r){
    if(l==r) tr[u]={l,r,a[l],0};
    else{
        tr[u]={l,r};
        int mid=l+r>>1;
        build(u<<1,l,mid);
        build(u<<1|1,mid+1,r);
        pushup(u);
    }
}

void modify(int u,int l,int r,int d){
    if(tr[u].l>=l and tr[u].r<=r) {
        tr[u].add+=d;
        tr[u].sum+=(LL)(tr[u].r-tr[u].l+1)*d;
    }
    else{
        pushdown(u);
        int mid=tr[u].l+tr[u].r>>1;
        if(l<=mid) modify(u<<1,l,r,d);
        if(r>mid) modify(u<<1|1,l,r,d);
        pushup(u);
    }
}

Node query(int u,int l,int r){
    if(tr[u].l>=l and tr[u].r<=r) return tr[u];
    else{
        pushdown(u);
        int mid=tr[u].l+tr[u].r>>1;
        if(r<=mid) return query(u<<1,l,r);
        else if(l>mid) return query(u<<1|1,l,r);
        else{
            Node left=query(u<<1,l,r);
            Node right=query(u<<1|1,l,r);
            Node res;
            pushup(res,left,right);
            return res;
        }
    }
}

int main()
{
    cin>>n>>m;
    for(int i=1;i<=n;i++) cin>>a[i];
    build(1,1,n);

    char op[2];
    int l,r,d;
    while(m--)
    {
        cin>>op>>l>>r;
        if(op[0]=='Q')
        {
            cout<<query(1,l,r).sum<<endl;
        }
        else
        {
            cin>>d;
            modify(1,l,r,d);
        }
    }
}

//fhq-treap(平衡树)

#include <bits/stdc++.h>
#define fastio ios::sync_with_stdio(0);cin.tie(0)
#define endl '\n'
using namespace std;
const int N=100010;
int root,idx;
int x,y,z;
struct Node
{
    int l,r;
    int key,val;
    int size;
}tr[N];

void pushup(int u)
{
    tr[u].size=tr[tr[u].l].size+tr[tr[u].r].size+1;
}

int getnode(int c)
{
    tr[++idx].key=c;
    tr[idx].val=rand();
    tr[idx].size=1;
    return idx;
}

void split(int u,int key,int &x,int &y)
{
    if(!u)
    {
        x=y=0;
        return;
    }
    else
    {
        if(tr[u].key<=key) x=u,split(tr[u].r,key,tr[u].r,y);
        else y=u,split(tr[u].l,key,x,tr[u].l);
    }
    pushup(u);
}

int merge(int a,int b)
{
    if(a==0 and b==0) return 0;
    if(tr[a].val > tr[b].val)
    {
        tr[a].r=merge(tr[a].r,b);
        pushup(a);
        return a;
    }
    else 
    {
        tr[b].l=merge(a,tr[b].l);
        pushup(b);
        return b;
    }
}

void insert(int c)
{
    split(root,c-1,x,y);
    root=merge(merge(x,getnode(c)),y);
}

void dele(int c)
{
    split(root,c,x,z);
    split(x,c-1,x,y);
    if(y) y=merge(tr[y].l,tr[y].r);
    root=merge(merge(x,y),z);
}

int getrank(int c)
{
    split(root,c-1,x,y);
    int ans=tr[x].size;
    root=merge(x,y);
    return ans;
}

int getnum(int c)
{
    int p=root;
    while(1)
    {
        if(tr[tr[p].l].size+1==c) break;
        else if(tr[tr[p].l].size+1 > c) p=tr[p].l;
        else c=c-tr[tr[p].l].size-1,p=tr[p].r;
    }
    return tr[p].key;
}

int getlast(int c)
{
    split(root,c-1,x,y);
    int p=x;
    while(tr[p].r) p=tr[p].r;
    int ans=tr[p].key;
    root=merge(x,y);
    return ans;
}

int getnext(int c)
{
    split(root,c,x,y);
    int p=y;
    while(tr[p].l) p=tr[p].l;
    int ans=tr[p].key;
    root=merge(x,y);
    return ans;
}

signed main()
{
    fastio;
    int n;
    cin>>n;
    for(int i=0;i<n;i++)
    {
        int k,x;
        cin>>k>>x;
        if(k==1) insert(x);
        else if(k==2) dele(x);
        else if(k==3) cout<<getrank(x)+1<<endl;
        else if(k==4) cout<<getnum(x)<<endl;
        else if(k==5) cout<<getlast(x)<<endl;
        else cout<<getnext(x)<<endl;
    }
}

//基础莫队

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 50010, M = 200010, S = 1000010;

int w[N];
int n,m,len;
int cnt[S];
int ans[M];
int res;

struct Query
{
    int id;
    int l;
    int r;
}q[M];

int get(int x)
{
    return x/len;
}

bool cmp(Query a,Query b)
{
    int i=get(a.l),j=get(b.l);
    if(i!=j) return i<j;
    return a.r<b.r;
}

void add(int x)
{
    if(!cnt[x]) res++;
    cnt[x]++;
}

void del(int x)
{
    cnt[x]--;
    if(!cnt[x]) res--;
}

int main()
{
    cin>>n;
    for(int i=1;i<=n;i++) cin>>w[i];
    
    cin>>m;
    
    len = max(1, (int)sqrt((double)n * n / m));
    
    for(int i=0;i<m;i++)
    {
        int l,r;
        cin>>l>>r;
        q[i]={i,l,r};
    }
    
    sort(q,q+m,cmp);
    
    for (int k = 0, i = 0, j = 1; k < m; k ++ )
    {
        int id = q[k].id, l = q[k].l, r = q[k].r;
        while (i < r) add(w[ ++ i]);
        while (i > r) del(w[i -- ]);
        while (j < l) del(w[j ++ ]);
        while (j > l) add(w[ -- j]);
        ans[id] = res;
    }

    
    for(int i=0;i<m;i++) cout<<ans[i]<<endl;
}







