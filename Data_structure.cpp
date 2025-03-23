//ST表

struct ST{
	std::vector<std::vector<int>> f;
    int n=0,m=0,flag=0;
    ST(){};
	ST(int nn,std::vector<int> &a,bool t=1){
        flag=t,n=nn,m=std::__lg(nn)+1;
        f.resize(n+1,std::vector<int> (m+1));
		for(int len=0;len<m;len++){
			for(int i=1;i+(1<<len)-1<=n;i++){
				if(!len) f[i][len]=a[i];
				else {
                    if(flag) f[i][len]=std::max(f[i][len-1],f[i+(1<<len-1)][len-1]);
                    else f[i][len]=std::min(f[i][len-1],f[i+(1<<len-1)][len-1]);
                }
			}
		}
	}
	int query(int l,int r){
		int k=std::__lg(r-l+1);
		if(flag) return std::max(f[l][k],f[r-(1<<k)+1][k]);
        return std::min(f[l][k],f[r-(1<<k)+1][k]);
	}
};

//并查集

struct DSU{
    std::vector<int> p, sz,add;
    DSU(int n): p(n), sz(n, 1),add(n,0){
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x){
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
        if(size(x) < size(y)) std::swap(x,y);
        sz[x] += sz[y];
        p[y] = x;
        add[y]-=add[x];
        return 1;
    }
    int size(int x){
        return sz[find(x)];
    }
};

struct DSU_Weight{//带权DSU
	std::vector<int> p, dist;
	DSU_Weight(int n): p(n), dist(n,0){
		iota(p.begin(), p.end(), 0);
	}
	int find(int x){
		if(p[x]!=x){
			int u=p[x];
			p[x]=find(p[x]);
			dist[x]+=dist[u];
		}
		return p[x];
	}
	bool merge(int x, int y,int w){
		int fx = find(x);
		int fy = find(y);
		if (fx == fy) return dist[x]==dist[y]+w;
		p[fx] = fy;
		dist[fx]=w+dist[y]-dist[x];
		return 1;
	}
};

//可撤销并查集

struct DSU{
    std::vector<int> p, sz;
    DSU(int n): p(n), sz(n, 1){
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x){
        return x == p[x] ? x : find(p[x]);
    }
    std::vector<std::array<int,3>> tmp;

    bool merge(int x, int y){
        x = find(x);
        y = find(y);
        if (x == y) return 0;
        if(size(x) < size(y)) std::swap(x,y);
        tmp.pb({y,x,sz[x]});
        sz[x] += sz[y];
        p[y] = x;
        return 1;
    }
    int size(int x){
        return sz[find(x)];
    }
    void roll(){
        while(tmp.size()){
            auto [y,x,t]=tmp.back();
            tmp.pop_back();
            p[y]=y;
            sz[x]=t;
        }
    }
};

//树状数组
struct BIT{
    int mn;
    std::vector<int> a;
    BIT(){};
    BIT(int n){
        a.clear(),a.resize(n+2);
        mn=n;
    };
    int lowbit(int x){return x&-x;}
    void add(int x,int c) {for(int i=x;i<=mn;i+=lowbit(i)) a[i]+=c;}
    long long sum(int x){
        long long res=0;
        for(int i=x;i;i-=lowbit(i)) res+=a[i];
        return res;
    }
};

//动态二维前缀和单点查询

struct presum_dynamic{
    int n;
	std::vector<std::vector<int>> a;
	void init(int nn){
        n=nn;a.clear();
		a.resize(n+2,vector<int> (n+2,0));
	}
	void updata(int x,int y,int d){
		for(int i=x;i<=n;i+=(i&-i)){
			for(int j=y;j<=n;j+=(j&-j)){
				a[i][j]+=d;
			}
		}
	}
	long long query(int x, int y){
		long long res = 0;
		for(int i=x;i;i-=(i&-i)){
			for(int j=y;j;j-=(j&-j)){
				res+=(long long)a[i][j];
			}
		}
		return res;
	}
	void add(int a,int b,int c,int d,int k){
		updata(a,b,k);
		updata(a,d+1,-k);
		updata(c+1,b,-k);
		updata(c+1,d+1,k);
	}
}tr;

//动态二维前缀和(区间修改区间查询)

template<class T>
struct presum_dynamic{
    std::vector<std::vector<T>> c1, c2, c3, c4;
    int szx,szy;
    presum_dynamic() {}
    presum_dynamic(int n,int m){init(n, m);}
    void init(int n, int m) {
        szx=n;szy=m;
    	c1.clear(),c2.clear(),c3.clear(),c4.clear();
        c1.resize(n+2,std::vector<T> (m+2,0));
        c2.resize(n+2,std::vector<T> (m+2,0));
        c3.resize(n+2,std::vector<T> (m+2,0));
        c4.resize(n+2,std::vector<T> (m+2,0));
    }
    void modify(int x, int y, T s) {
        for(int i=x;i<=szx;i+=(i&-i)){
            for(int j=y;j<=szy;j+=(j&-j)) {
                c1[i][j]+=s;
                c2[i][j]+=s*y;
                c3[i][j]+=s*x;
                c4[i][j]+=s*x*y;
            }
        }
    }
    T query(int x, int y) {
        T s = 0;
        for(int i=x;i;i-=(i&-i)){
            for(int j=y;j;j-=(j&-j)) {
                s+=c1[i][j]*(x+1)*(y+1);
                s-=c2[i][j]*(x+1);
                s-=c3[i][j]*(y+1);
                s+=c4[i][j];
            }
        }
        return s;
    }
    void recmodify(int a, int b, int c, int d,T s) {
        modify(a,b,s);
        modify(a,d+1,-s);
        modify(c+1,b,-s);
        modify(c+1,d+1,s);
    }
    T recquery(int a, int b, int c, int d) {
        return query(c,d)-query(a-1,d)-query(c,b-1)+query(a-1,b-1);
    }
};

//笛卡尔树(max/min,子树连续)

struct Descar_tree{
	std::stack<int> stk;
	std::vector<std::array<int,2>> tr;
    std::vector<int> L,R;
	int root=-1;
    Descar_tree(){};
    Descar_tree(int n,std::vector<int> &a){
        L.clear(),R.clear();
        L.resize(n+2),R.resize(n+2);
        tr.resize(n+2);root=-1;
        stk=std::stack<int>();
		for(int i=1;i<=n;i++){
			while(stk.size() and a[stk.top()]<a[i]) {
				tr[i][0]=stk.top(),stk.pop();
			}
			if(stk.size()) tr[stk.top()][1]=i;
			else root=i;
			stk.push(i);
		}
        dfs(root,-1);
    };
    void dfs(int u,int fa){
        auto &t=tr[u];
        L[u]=R[u]=u;
        if(t[0]){
            dfs(t[0],u);
            L[u]=std::min(L[u],L[t[0]]);
            R[u]=std::max(R[u],R[t[0]]);
        }
        if(t[1]){
            dfs(t[1],u);
            L[u]=std::min(L[u],L[t[1]]);
            R[u]=std::max(R[u],R[t[1]]);
        }
    };
};

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

//线段树：//根据题目修改build函数与Tag、Info对应apply函数

struct Tag{//lazy tag
	Mint add=0,mul=1;
	Tag(Mint a=0,Mint b=1):add(a),mul(b) {};
	void apply(Tag t){
		add=(add*t.mul+t.add);
		mul=(mul*t.mul);
	};
};

struct Info{//information
	int l,r;
	Mint sum;
	Info() {};
	void apply(Tag t){
		sum=sum*t.mul+(r-l+1)*t.add;
	};
	friend Info operator + (const Info &a,const Info &b){
		Info res;
		res.l=a.l,res.r=b.r;
		res.sum=a.sum+b.sum;
		return res;
	};
};

struct Segment_tree{
	Info info[N*4];
	Tag tag[N*4];
	
	void pushup(int u){
		info[u]=info[u<<1]+info[u<<1|1];
	}
	
	void apply(int u,const Tag &t){
		info[u].apply(t);
		tag[u].apply(t);
	}
	
	void pushdown(int u){
		apply(u<<1,tag[u]);
		apply(u<<1|1,tag[u]);
		tag[u]=Tag();
	}
	
	void build(int u,int l,int r){
		tag[u]=Tag();
		info[u].l=l,info[u].r=r;
		if(l==r) {
			info[u].l=l,info[u].r=r;
			info[u].sum=a[l];
		}else{
			int mid=(l+r)>>1;
			build(u<<1,l,mid);
			build(u<<1|1,mid+1,r);
			pushup(u);
		}
	}
	
	void modify(int u,int l,int r,int pl,int pr,const Tag &t){
		if(l>=pl and r<=pr) {
			apply(u,t);
			return;
		}
		pushdown(u);
		int mid=(l+r)>>1;
		if(pl<=mid) modify(u<<1,l,mid,pl,pr,t);
		if(pr>mid) modify(u<<1|1,mid+1,r,pl,pr,t);
		pushup(u);
	}
	
	Info query(int u,int l,int r,int pl,int pr){
		if(l>=pl and r<=pr) return info[u];
		pushdown(u);
		int mid=(l+r)>>1;
		if(pl>mid) return query(u<<1|1,mid+1,r,pl,pr);
		else if(pr<=mid) return query(u<<1,l,mid,pl,pr);
		else return query(u<<1,l,mid,pl,pr)+query(u<<1|1,mid+1,r,pl,pr);
	}
}tr;

//动态开点版本线段树

struct Segment_tree{
	std::vector<Info> info;
	std::vector<Tag> tag;
	std::vector<int> ls,rs;
	int idx=-1,root;
	
	int getnode(int l,int r,int v){
		idx++;
		info.push_back(Info(v));
		ls.push_back(-1),rs.push_back(-1);
		tag.push_back(Tag());
		return idx;
	}
	
	void pushup(int u){
		info[u]=info[ls[u]]+info[rs[u]];
	}
	
	void apply(int u,int l,int r,const Tag &t){
		info[u].apply(t,l,r);
		tag[u].apply(t);
	}
	
	void pushdown(int u,int l,int r){
		if(ls[u]==-1){
			int t=getnode(l,(l+r)/2,0);
			ls[u]=t;
		}
		if(rs[u]==-1){
			int t=getnode((l+r)/2+1,r,0);
			rs[u]=t;
		}
		apply(ls[u],l,(l+r)/2,tag[u]);
		apply(rs[u],(l+r)/2+1,r,tag[u]);
		tag[u]=Tag();
	}
	
	void modify(int u,int l,int r,int pl,int pr,const Tag &t){
		if(l>=pl and r<=pr) {
			apply(u,l,r,t);
			return;
		}
		pushdown(u,l,r);
		int mid=(l+r)>>1;
		if(pl<=mid) modify(ls[u],l,mid,pl,pr,t);
		if(pr>mid) modify(rs[u],mid+1,r,pl,pr,t);
		pushup(u);
	}
	
	void init(){
		root=getnode(1,n,0);
		for(int i=1;i<=n;i++) modify(0,1,n,i,i,Tag(a[i]));
	}
	
	Info query(int u,int l,int r,int pl,int pr){
		if(l>=pl and r<=pr) return info[u];
		pushdown(u,l,r);
		int mid=(l+r)>>1;
		if(pl>mid) return query(rs[u],mid+1,r,pl,pr);
		else if(pr<=mid) return query(ls[u],l,mid,pl,pr);
		else return query(ls[u],l,mid,pl,pr)+query(rs[u],mid+1,r,pl,pr);
	}
}tr;

//可持久化线段树(Persistent Segment Tree)求区间小于等于k数量

int root[N],idx;
struct Node{
	int ls,rs;
	int v;
}tr[N*20];

int build(int l,int r){
	int u=++idx;
	tr[u]={-1,-1,0};
	if(l==r) return u;
	int mid=(l+r)/2;
	tr[u].ls=build(l,mid);
	tr[u].rs=build(mid+1,r);
	return u;
}

int update(int p,int x,int c,int l,int r){
	int u=++idx;
	tr[u]=tr[p];
	tr[u].v+=c;
	if(l==r) return u;
	int mid=(l+r)/2;
	if(x<=mid) tr[u].ls=update(tr[u].ls,x,c,l,mid);
	else tr[u].rs=update(tr[u].rs,x,c,mid+1,r);
	return u;
}

int query(int p,int q,int k,int l,int r){
	if(l==r) return tr[q].v-tr[p].v;
	int mid=(l+r)/2;
	if(k<=mid) return query(tr[p].ls,tr[q].ls,k,l,mid);
	else {
		int sum=tr[tr[q].ls].v-tr[tr[p].ls].v;
		return sum+query(tr[p].rs,tr[q].rs,k,mid+1,r);
	}
	return 0;
}

//fhq-treap(平衡树)

struct FHQ_treap{
	#define Maxn 100010
	#define ls tree[p].pl
	#define rs tree[p].pr
	private:
		int All=0,root=0;
		struct NODE { int pl,pr,siz,cnt,rnd; int val; };
		NODE tree[Maxn];
		inline int Dot() { return ++All; }
		inline int New(int Val){
			int p=Dot();
			tree[p].rnd=rand(),tree[p].val=Val;
			tree[p].siz=tree[p].cnt=1;
			tree[p].pl=tree[p].pr=0;
			return p;
		}
		inline void pushdown(int p) { p--; }
		inline void pushup(int p) { 
		    tree[p].siz=tree[ls].siz+tree[rs].siz+tree[p].cnt; 
		}
		void split(int p,int k,int &x,int &y){
			if(!p) { x=y=0; return; }
			pushdown(p);
			if(tree[p].val<=k) x=p,split(rs,k,rs,y);
			else y=p,split(ls,k,x,ls);
			pushup(p);
		}
		int merge(int x,int y){
			if(!x || !y) return x+y;
			if(tree[x].rnd<tree[y].rnd){
				pushdown(x),tree[x].pr=merge(tree[x].pr,y),pushup(x);
				return x;
			}else{
				pushdown(y),tree[y].pl=merge(x,tree[y].pl),pushup(y);
				return y;
			}
		}
		inline int kth(int p,int Rank){
			while(p){
				if(tree[ls].siz>=Rank) p=ls;
				else if(tree[ls].siz+tree[p].cnt>=Rank) return p;
				else Rank-=tree[ls].siz+tree[p].cnt,p=rs;
			}
			return p;
		}
		int x,y,z;
	public:
		void Insert(int Val){
		    split(root,Val,x,y),root=merge(merge(x,New(Val)),y); 
		}
		void Delete_one(int Val){
			split(root,Val,x,z),split(x,Val-1,x,y);
			y=merge(tree[y].pl,tree[y].pr);
			root=merge(merge(x,y),z);
		}
		int Rank_to_Value(int Rank)
			{ return tree[kth(root,Rank)].val; }
		int Value_to_Rank(int Value){
			split(root,Value-1,x,y);
			int ret=tree[x].siz+1;
			root=merge(x,y);
			return ret;
		}
		int Findpre(int Value){
			split(root,Value-1,x,y);
			int ret=tree[kth(x,tree[x].siz)].val;
			root=merge(x,y);
			return ret;
		}
		int Findnex(int Value){
			split(root,Value,x,y);
			int ret=tree[kth(y,1)].val;
			root=merge(x,y);
			return ret;
		}
}tr;


//莫队

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;


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







