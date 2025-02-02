//最短路dijkstra

std::vector<int> dijkstra(int S){
	std::vector<int> dist(n+1,-1);
	std::vector<bool> st(n+1,0);
	std::priority_queue<std::array<int,2>,std::vector<std::array<int,2>>,std::greater<std::array<int,2>>> q;
	q.push({0,S});dist[S]=0;
	while(q.size()){
		auto t=q.top();
		q.pop();
		if(st[t[1]]) continue;
		st[t[1]]=1;
		for(auto x:eg[t[1]]){
			if(dist[x]==-1 or dist[x]>t[0]+1){
				dist[x]=t[0]+1;
				q.push({dist[x],x});
			}
		}
	}
	return dist;
}


//可重复走边非严格次短路dijkstra

void dijkstra2(int u){
    for (int i=0; i<N;i++) dis1[i]=dis2[i]=1e18;
    std::priority_queue<std::array<int,2>,std::vector<std::array<int,2>>,std::greater<>> q;
    q.push({0,u});
    dis1[u]=0;
    while (!q.empty()){
        auto t = q.top();
        q.pop();
        if (dis2[t[1]]<t[0]) continue;
        for (auto [j,w]:eg[t[1]]){
            int nw = t[0] + w;
            if(nw<dis1[j]){
                std::swap(dis1[j],nw);
                q.push({dis1[j], j});
            }
            if(nw<dis2[j]) {
                std::swap(dis2[j],nw);
                q.push({dis2[j],j});
            }
        }
    }
}

//可重复走边严格次短路

void dijkstra2(int u){
    for (int i=0; i<N;i++) dis1[i]=dis2[i]=1e18;
    std::priority_queue<std::array<int,2>,std::vector<std::array<int,2>>,std::greater<>> q;
    q.push({0,u});
    dis1[u]=0;
    while (!q.empty()){
        auto t = q.top();
        q.pop();
        if (dis2[t[1]]<t[0]) continue;
        for (auto [j,w]:eg[t[1]]){
            int nw = t[0] + w;
            if(nw<dis1[j]){
                std::swap(dis1[j],nw);
                q.push({dis1[j], j});
            }
            if(nw<dis2[j] and dis1[j]!=nw) {
                std::swap(dis2[j],nw);
                q.push({dis2[j],j});
            }
        }
    }
}

//点分治处理树上路径相关计数

bool st[N];
int get_size(int u,int f){
	if(st[u]) return 0;
	int sz=1;
	for(auto &[x,w]:eg[u]){
		if(x==f) continue;
		sz+=get_size(x,u);
	}
	return sz;
}

int getwc(int u,int f,int cnt,int &root){
	if(st[u]) return 0;
	int sz=1,ma=0;
	for(auto &[c,w]:eg[u]){
		if(c==f) continue;
		int t=getwc(c,u,cnt,root);
		sz+=t;
		ma=max(t,ma);
	}
	ma=max(ma,cnt-sz);
	if(ma<=cnt/2) root=u;
	return sz;
}

void dfs(int u,int f){
	if(st[u]) return;
	getwc(u,0,get_size(u,0),u);
	st[u]=1;
	int idx=0;

    int cnt[2]={0};

	for(auto &[x,w]:eg[u]){
		if(st[x] or x==f) continue;
        int q[2]={0};
		function<void(int,int,int)> get_dist=[&](int uu,int ff,int dis){
			if(st[uu]) return;
            q[dis]++;
			for(auto &[xx,ww]:eg[uu]){
				if(xx==ff) continue;
				get_dist(xx,uu,(dis+ww)&1);
			}
		};
		get_dist(x,u,w&1);
        res+=q[1];
        res+=q[1]*cnt[0]+cnt[1]*q[0];
        cnt[1]+=q[1],cnt[0]+=q[0];
	}
	for(auto &[x,w]:eg[u]){
		if(x==f) continue;
		dfs(x,u);
	}
}

//kruskal最小生成树

struct kurskal{
	std::vector<std::vector<std::array<int,2>>> g;
	int cnt,sum;
	
	void build(int n,std::vector<std::array<int,3>> edge){
		g.clear();g.resize(n+1);
		cnt=0,sum=0;
		sort(edge.begin(),edge.end(),[&](auto l,auto r){
			return l[2]<r[2];
		});
		DSU dsu(n);
		for(auto [u,v,w]:edge){
			if(dsu.same(u,v)) continue;
			sum+=w;
			dsu.merge(u,v);
			g[u].pb({v,w});
			g[v].pb({u,w});
		}
		cnt=dsu.size(1);
	}
}kt;

//kruskal重构树

int tr[N*2],root;
void kruskal_tree(){
	sort(eg,eg+cnt);
	for(int i=0;i<N*2;i++) p[i]=i;
	memset(h,-1,sizeof h);idx=0;
	for(int i=0;i<cnt;i++){
		int a=find(eg[i].a),b=find(eg[i].b);
		if(a==b) continue;
		root++;
		tr[root]=eg[i].h;
		merge(a,root);
		merge(b,root);
		add(root,a,0);add(root,b,0);
	}
}

//倍增LCA（最近公共祖先）

struct LCA{
    std::vector<int> dep;
    std::vector<std::vector<int>> eg,fa;
    int m;
    void init(int n){
        m=log(n)/log(2)+1;
        eg.clear(),dep.clear(),fa.clear();
        eg.resize(n+1),dep.resize(n+1,0);
        fa.resize(n+1,vector<int> (m+1,0));
    }
    
    void dfs(int u, int f) {
        for (int i = 1; i <= m; i++) {
            fa[u][i] = fa[fa[u][i - 1]][i - 1];
        }
        for (auto j:eg[u]) {
            if (j == f) continue;
            fa[j][0] = u;
            dep[j] = dep[u] + 1;
            dfs(j, u);
        }
    }
    
    int lca(int a, int b) {
        if (dep[a] < dep[b]) swap(a, b);
        for (int i = m; i >= 0; i--) {
            if (dep[a] - (1 << i) >= dep[b]) a = fa[a][i];
        }
        if (a == b) return a;
        for (int i = m; i >= 0; i--) {
            if (fa[a][i] != fa[b][i]) {
                a = fa[a][i];
                b = fa[b][i];
            }
        }
        return fa[a][0];
    }
}tr;

//树链剖分，修改路径，修改子树

struct HLD{
    int n,cnt;
    std::vector<int> sz,dep,top,son,fa,in,out,seq;
    std::vector<std::vector<int>> g;
    void init(int nn,std::vector<std::vector<int>> &eg,int root){
        n=nn;cnt=0;g=eg;
        sz.clear(),dep.clear();
        top.clear(),seq.clear();
        son.clear(),fa.clear();
        in.clear(),out.clear();
        sz.resize(n+1),dep.resize(n+1);
        top.resize(n+1),seq.resize(n+1);
        son.resize(n+1),fa.resize(n+1);
        in.resize(n+1),out.resize(n+1);
        
        dep[root] = 1;
        top[root] = root;
        dfs1(root);
        dfs2(root);
    };
    void dfs1(int u){
        if (fa[u]) g[u].erase(find(g[u].begin(),g[u].end(),fa[u]));
        sz[u]=1;
        for(auto &j:g[u]){
            fa[j]=u;
            dep[j]=dep[u]+1;
            dfs1(j);
            sz[u]+=sz[j];
            if (sz[j]>sz[g[u][0]])
                std::swap(j,g[u][0]);
        }
    }
    void dfs2(int u){
        in[u]=++cnt;
        seq[in[u]]=u;
        a[in[u]]=w[u];
        //do something
        for (auto j : g[u]){
            top[j]=(j == g[u][0] ? top[u] : j);
            dfs2(j);
        }
        out[u]=cnt;
    }
    int lca(int u,int v){
        while(top[u]!=top[v]){
            if(dep[top[u]]<dep[top[v]]) v=fa[top[v]];
            else u=fa[top[u]];
        }
        return dep[u]<dep[v] ? u : v;
    }
    int dist(int u,int v){
        return dep[u]+dep[v]-2*dep[lca(u,v)];
    }
    std::vector<std::array<int,2>> get_path(int u,int v,int is_edge=0){
        std::vector<std::array<int,2>> res;
        while(top[u]!=top[v]){
            if(dep[top[u]]<dep[top[v]]) std::swap(u,v);
            res.push_back({in[top[u]],in[u]});
            u=fa[top[u]];
        }
        if(dep[u]>dep[v]) std::swap(u,v);
        if(is_edge and u==v) return res;
        res.push_back({in[u]+is_edge,in[v]});
        return res;
    }
    std::vector<std::array<int,2>> get_tree(int u,int is_edge=0){
    	if(is_edge and in[u]==out[u]) assert(0);
        return {{in[u]+is_edge,out[u]}};
    }
}hld;

//二分图最大匹配HK

struct HopcroftKarp {//idx_start_1
	std::vector<int> g,l,r;
	int ans;
	HopcroftKarp(int n,int m,const std::vector<std::array<int,2>> &e): g(e.size()+2),l(n+2, -1),r(m+2, -1),ans(0){
		std::vector<int> deg(n+2);
		for(auto &[x, y]:e) deg[x]++;
		for(int i=1;i<=n+1;i++) deg[i]+=deg[i-1];
		for(auto &[x, y]:e) g[--deg[x]]=y;
		std::vector<int> q(n+2);
		while(true){
			std::vector<int> a(n+2,-1),p(n+2,-1);
			int t=0;
			for(int i=1;i<=n;i++) if(l[i]==-1)
				q[t++]=a[i]=p[i]=i;
			bool match=false;
			for(int i=0;i<t;i++) {
				int x=q[i];
				if(~l[a[x]]) continue;
				for(int j=deg[x];j<deg[x+1];j++) {
					int y=g[j];
					if(r[y]==-1) {
						while(~y) r[y]=x,std::swap(l[x],y),x=p[x];
						match=true,ans++;
						break;
					}
					if(p[r[y]]==-1)
						q[t++]=y=r[y],p[y]=x,a[y]=a[x];
				}	
			}
			if(!match) break;
		}
	}
};

//最大流:节点编号(1-n)

struct MF {
  int h[N], e[M], ne[M], w[M], idx = 0;
  int vis[N];
  
  int n, S, T;
  int maxflow = 0;
  int dep[N], cur[N];

  void init(int n_,int S_,int T_) {
    n=n_,S=S_,T=T_;
    maxflow=0;
    memset(h, -1, sizeof h);
    idx = 0;
  }

  void add(int a, int b, int c) {
    e[idx] = b, ne[idx] = h[a], w[idx] = c, h[a] = idx++;
    e[idx] = a, ne[idx] = h[b], w[idx] = 0, h[b] = idx++;
  }

  bool bfs() {
    std::queue<int> q;
    memset(dep, -1, sizeof dep);
    dep[S] = 0;
    q.push(S);
    cur[S] = h[S];
    while (q.size()) {
      int u = q.front();
      q.pop();
      for (int i = h[u]; ~i; i = ne[i]) {
        int j = e[i];
        if (dep[j] == -1 and w[i]) {
          dep[j] = dep[u] + 1;
          cur[j] = h[j];
          if (j == T) return 1;
          q.push(j);
        }
      }
    }
    return 0;
  }

  int dfs(int u, int limit) {
    if (u == T or !limit) return limit;
    int flow = 0;
    for (int i = cur[u]; ~i and flow < limit; i = ne[i]) {
      cur[u] = i;
      int j = e[i];
      if (dep[j] == dep[u] + 1 and w[i]) {
        int t = dfs(j, std::min(limit - flow, w[i]));
        if (!t) dep[j] = -1;
        w[i] -= t;
        w[i ^ 1] += t;
        flow += t;
      }
    }
    return flow;
  }

  void dinic() {
    while (bfs()) {
      memcpy(cur, h, sizeof h);
      maxflow += dfs(S, INF);
    }
  }
  void dfs(int u) {//找与源点同一侧的点
	  vis[u] = 1;
	  for (int i = h[u]; ~i; i = ne[i]) {
	    int v = e[i];
	    if (!vis[v] && w[i]) dfs(v);
	  }
  }
} mf;

//最小费用最大流:节点编号(1-n)

namespace dinic{
	const int N=1e5+7,M=2e6+7;
	const int INF=1e9;
	int n,S,T;
	int head[N],ver[M],nex[M],tot,cur[N];
	int dist[N],edge[M],cost[M],maxflow,mincost;
	bool vis[N];
	
	inline void add(int x,int y,int z,int c=0,bool o=1){
		ver[tot]=y;
		edge[tot]=z;
		cost[tot]=c;
		nex[tot]=head[x];
		head[x]=tot++;
		if(o) add(y,x,0,-c,0);
	}
	inline bool spfa(){
		for(int i=1;i<=n;i++) dist[i]=INF;
		memset(vis,0,sizeof vis);
		std::queue<int> q;
		q.push(S);
		dist[S]=0,vis[S]=1;
		while(q.size()){
			auto x=q.front();
			q.pop();
			vis[x]=0;
			for(int i=head[x];~i;i=nex[i]){
				int y=ver[i];
				int z=edge[i],c=cost[i];
				if(dist[y]<=dist[x]+c or !z) continue;
				dist[y]=dist[x]+c;
				if(!vis[y]) q.push(y),vis[y]=1;
			}
		}
		return dist[T]!=INF;
	}
	int dfs(int x,int flow=INF){
		if(x==T) return flow;
		int ans=0,k,i;
		vis[x]=1;
		for(int i=cur[x];~i and flow;i=nex[i]){
			cur[x]=i;
			int y=ver[i];
			int z=edge[i],c=cost[i];
			if(!z or (dist[y]!=dist[x]+c) or vis[y]) continue;
			k=dfs(y,std::min(flow,z));
			if(!k) dist[y]=INF;
			edge[i]-=k;
			edge[i^1]+=k;
			ans+=k;
			mincost+=k*c;
			flow-=k;
		}
		vis[x]=0;
		return ans;
	}
	inline void main(){
		while(spfa()){
			for(int i=1;i<=n;i++){
				cur[i]=head[i];
			}
			int now;
			while((now=dfs(S,INF))) maxflow+=now;
		}
	}
	inline void init(int _n,int _S,int _T){
		n=_n,S=_S,T=_T;
		tot=0,maxflow=0,mincost=0;
		memset(head,-1,sizeof head);
	}
}

// Tarjan求边双连通分量E_DCC

struct E_DCC{
	int n,cnt,cnt_edge;
	int dfn[N], low[N];
	std::vector<std::pair<int, int>> e[N];
	std::vector<std::vector<int>> ans;
	std::stack<int> st;
	std::vector<std::array<int,2>> bridge;
    int is_bridge[N];
	
	void init(int nn){
		n=nn,cnt_edge=cnt=0;
		st=std::stack<int>();
		ans.clear();bridge.clear();
		for(int i=1;i<=n;i++) {
			e[i].clear();
			dfn[i]=low[i]=0;
			is_bridge[i]=0;
		}
	}
	
	void add(int u,int v){
		cnt_edge++;
		e[u].push_back(std::make_pair(v, cnt_edge<<1));
		e[v].push_back(std::make_pair(u, cnt_edge<<1|1));
	}
	
	void dfs(int x, int las){
		low[x] = dfn[x] = ++cnt;
		st.push(x);
		for (auto i: e[x]){
			if (i.second == (las ^ 1)) continue;
			if (!dfn[i.first]){
				dfs(i.first, i.second);
				low[x] = std::min(low[x], low[i.first]);
				if(low[i.first]>dfn[x]){
                    bridge.push_back({x,i.first});
                    is_bridge[x]=1;
                    is_bridge[i.first]=1;
                }
			}else low[x] = std::min(low[x], dfn[i.first]);
		}
		if (dfn[x] == low[x]){
			std::vector<int> vec;
			vec.push_back(x);
			while (st.top() != x){
				vec.push_back(st.top());
				st.pop();
			}
			st.pop();
			ans.push_back(vec);
		}
	}
	
	void tarjan(){
		for(int i=1;i<=n;i++) 
			if(!dfn[i]) dfs(i,0);
	}
}tj;

// Tarjan求点双连通分量V_DCC

struct V_DCC{
    std::vector<int> eg[N];
    int dfn[N],low[N],timespace;
    int stk[N],top;
    bool cut[N];
    int dcc_cnt,n;
    std::vector<int> dcc[N];
    void init(int x){
        timespace=1;
        n=x;
    }
    void add(int a,int b){
        eg[a].push_back(b);
        eg[b].push_back(a);
    }
    void dfs(int u,int root){
        dfn[u]=low[u]=timespace++;
        stk[++top]=u;
        
        if(root==u and eg[u].empty()){
            dcc_cnt++;
            dcc[dcc_cnt].push_back(u);
            return;
        }
        int cnt=0;
        for(auto j:eg[u]){
            if(!dfn[j]){
                dfs(j,root);
                low[u]=std::min(low[u],low[j]);
                if(dfn[u] <= low[j]){
                    cnt++;
                    if(u!=root or cnt > 1) cut[u]=1;
                    dcc_cnt++;
                    int y;
                    do{
                        y=stk[top--];
                        dcc[dcc_cnt].push_back(y);
                    }while(y!=j);
                    dcc[dcc_cnt].push_back(u);
                }
            }
            else low[u]=std::min(low[u],dfn[j]);
        }
    }
    void tarjan(){
        for(int i=1;i<=n;i++){
            if(!dfn[i]) dfs(i,i);
        }
    }
}tj;

// Tarjan求有向图强连通分量SCC

struct Tarjan{
    std::vector<int> eg[N];
    int dfn[N],low[N],timespace;
    int stk[N],top;
    bool st[N];
    int color[N];
    int scc_cnt,n;
    std::vector<int> scc[N];
    void init(int x){
        timespace=1;scc_cnt=0;
        top=0;
        n=x;
        for(int i=1;i<=n;i++) eg[i].clear(),scc[i].clear();
        for(int i=1;i<=n;i++){
            dfn[i]=0,low[i]=0;
            color[i]=0;
        }
    }
    void add(int a,int b){
        eg[a].push_back(b);
    }
    void dfs(int u){
        low[u]=dfn[u]=++timespace;
        stk[++top]=u;
        st[u]=1;
        for(auto v:eg[u]){
            if (!dfn[v]){
                dfs(v);
                low[u]=std::min(low[u],low[v]);
            }else if(st[v]) low[u]=std::min(low[u],low[v]);
        }
        if(dfn[u]==low[u]){
            scc_cnt++;
            int y;
            while (y=stk[top--])
            {
                scc[scc_cnt].push_back(y);
                color[y]=scc_cnt;
                st[y]=0;
                if (u==y) break;
            }
        }
    }
    void tarjan(){
        for(int i=1;i<=n;i++){
            if(!dfn[i]) dfs(i);
        }
    }
}tj;



//线段树优化建图区间连边
struct segment_tree_graph{
	int cnt=0,cntin=0;
	int root_in,root_out;
	int in[N],out[N];
    array<int,2> tr[N*10];
	int h[N*10],e[N*20],ne[N*20],w[N*20],idx;
    int rd[N*10],st[N*10];
 
    bool is_init=0;
 
    void add_edge(int a,int b,int c){
        e[idx]=b,w[idx]=c,ne[idx]=h[a],h[a]=idx++;
        rd[b]++;
    }
    void build_in_tree(int u,int l,int r){
        if(l==r){
            in[l]=u;
            st[u]=1;
            cnt=max(cnt,u);
            return;
        }
        int mid=(l+r)/2;
        build_in_tree(u*2,l,mid);
        build_in_tree(u*2+1,mid+1,r);
        tr[u]={u*2,u*2+1};
        add_edge(u*2,u,0);
        add_edge(u*2+1,u,0);
    }
    void build_out_tree(int u,int l,int r){
        if(l==r){
            out[l]=u;
            return;
        }
        int mid=(l+r)/2;
        int ls=++cnt,rs=++cnt;
        build_out_tree(ls,l,mid);
        build_out_tree(rs,mid+1,r);
        tr[u]={ls,rs};
        add_edge(u,ls,0);
        add_edge(u,rs,0);
    }
	void build(int l,int r){
        if(!is_init){
            for(int i=0;i<N*10;i++) h[i]=-1;
            is_init=1;
        }
        root_in=++cnt;
		build_in_tree(root_in,l,r);
        cntin=cnt,root_out=++cnt;
        build_out_tree(root_out,l,r);
        for(int i=1;i<=n;i++) add_edge(out[i],in[i],0),add_edge(in[i],out[i],0);
	}
	
	void link(int u,int l,int r,int pl,int pr,int x,int w,int flag){
		if(l>=pl and r<=pr){
			if(flag) add_edge(x,u,0);
			else add_edge(u,x,w);
			return;
		}
		int mid=(l+r)/2;
		if(pl<=mid) link(tr[u][0],l,mid,pl,pr,x,w,flag);
		if(pr>mid) link(tr[u][1],mid+1,r,pl,pr,x,w,flag);
	}
	
	void modify(int l1,int r1,int l2,int r2,int w){
		int nt=++cnt;
        link(root_in,1,n,l1,r1,nt,w,0);
		link(root_out,1,n,l2,r2,nt,w,1);
		
	}
	void clear(){
		for(int i=0;i<=cnt;i++){
			tr[i]={0,0};
			in[i]=out[i]=0;
			h[i]=-1;
            rd[i]=0,st[i]=0;
		}
		idx=0,cnt=0;root_in=root_out=0;
	}
}tr;


//虚树dp

struct virtual_tree{
	vector<vector<int>> veg;
	vector<int> f,st;
	vector<int> stk;
	int top;
	void init(int n){
		veg.resize(n+1);
		f.resize(n+1,0);
		st.resize(n+1,0);
		top=0;
	}
	void clear(){
		stk.clear();top=0;
		function<void(int,int)> dclear=[&](int u,int la){
			for(auto x:veg[u]) if(x!=la) dclear(x,u);
			veg[u].clear();
			f[u]=0;st[u]=0;
		};
		dclear(1,0);
	}
	void add_edge(int a,int b){
		veg[a].push_back(b);
	}
	void build(vector<int> &q){
		sort(all(q),[&](int u,int v){
			return dfn[u]<dfn[v];
		});
		stk.push_back(0);
		stk.push_back(1),top++;
		for(auto x:q){
			st[x]=1;
			if(x==1) continue;
			int t=lca(x,stk[top]);
			if(t!=stk[top]){
				while(dfn[t]<dfn[stk[top-1]]){
					add_edge(stk[top-1],stk[top]);
					top--;stk.pop_back();
				}
				if(dfn[t]>dfn[stk[top-1]]){
					add_edge(t,stk[top]);
					stk.pop_back();
					stk.push_back(t);
				}else{
					add_edge(t,stk[top]);
					top--;stk.pop_back();
				}
			}
			top++;stk.push_back(x);
		}
		while(stk.size()>=2){
			auto t=stk.back();
			stk.pop_back();
			if(stk.back()==0) break;
			else add_edge(stk.back(),t);
		}
	}


	int dp(int u,int la){

	}
}vtr;

//欧拉回路

struct Euler_path{
	int S,T;
	bool check_graph(){
		vector<int> l,r;
		for(int i=1;i<=n;i++){
			if(rd[i]==cd[i]) continue;
			if(rd[i]==cd[i]+1) r.pb(i);
			else if(cd[i]==rd[i]+1) l.pb(i);
			else return 0;
		}
		if(l.empty() and r.empty()) S=1,T=2;
		else if(l.size()==1 and r.size()==1) S=l[0],T=r[0];
		else return 0;
		return 1;
	};
	
	vector<int> check_and_get_ugraph(vector<array<int,2>> &edge){
		vector<vector<array<int,2>>> g(n+2);
		vector<int> d(n+2,0),pt;
		int len=edge.size();
		for(int i=0;i<len;i++){
			auto [u,v]=edge[i];
			g[u].pb({v,i});
			g[v].pb({u,i});
			d[v]++,d[u]++;
		}
		
		for(int i=1;i<=n;i++) if(d[i]&1) pt.pb(i);
		if(pt.empty()) S=1;
		else if(pt.size()==2) S=pt[1];
		else return {};
		
		vector<int> st(len+2),res;
		function<void(int)> dfs=[&](int u){
			while(g[u].size()){
				auto [v,id]=g[u].back();
				g[u].pop_back();
				if(!st[id]){
					st[id]=1;
					dfs(v);
				}
			}
			res.pb(u);
		};
		dfs(S);
		reverse(all(res));
		return res;
	};
	
	vector<int> get_graph(){
		vector<int> res;
		function<void(int)> dfs=[&](int u){
			while(eg[u].size()){
				auto v=eg[u].back();
				eg[u].pop_back();
				dfs(v);
			}
			res.pb(u);
		};
		dfs(S);
		reverse(all(res));
		return res;
	};
}ep;

//支配树 推荐博客: https://zerol.me/2018/10/22/dominator-tree/

struct DominatorTree{
    int n,idx;
    std::vector<int> dfn,rev,fa,p,mn,idom,sdom;
    std::vector<std::vector<int>> eg,reg,tr;

    DominatorTree(){}
    DominatorTree(int n,std::vector<std::array<int,2>> &Edge,int root=-1): n(n) {
        assert(root!=-1);

        dfn.resize(n+2),rev.resize(n+2);
        fa.resize(n+2),p.resize(n+2);
        idom.resize(n+2),sdom.resize(n+2);
        eg.resize(n+2),reg.resize(n+2);
        mn.resize(n+2),tr.resize(n+2);

        idx=0;
        for(int i=0;i<=n;i++){
            p[i]=i,sdom[i]=i,mn[i]=i;
        }
        for(auto [u,v]:Edge) add(u,v);
        tarjan(root);
    };  
    void add(int a,int b){
        eg[a].pb(b);
        reg[b].pb(a);
    }
    void dfs(int u){
        dfn[u]=++idx;
        rev[idx]=u;
        for(auto x:eg[u]){
            if(!dfn[x]){
                fa[x]=u;
                dfs(x);
            }
        }
    }
    int find(int x){
        if(x==p[x]) return x;
        int j=find(p[x]);
        if(dfn[sdom[mn[p[x]]]]<dfn[sdom[mn[x]]]) mn[x]=mn[p[x]];
        return p[x]=j;
    }
    void tarjan(int root){
        dfs(root);
        for(int i=n;i>1;i--){
            int u=rev[i];
            for(auto x:reg[u]){
                find(x);
                if(dfn[sdom[mn[x]]]<dfn[sdom[u]]) sdom[u]=sdom[mn[x]];
            }
            p[u]=fa[u];
            tr[sdom[u]].pb(u);
            u=fa[u];
            for(auto x:tr[u]){
                find(x);
                idom[x]=(u==sdom[mn[x]]) ? u : mn[x];
            }
            tr[u].clear();
        }
        for(int i=2;i<=n;i++){
            if(idom[rev[i]]!=sdom[rev[i]]) idom[rev[i]]=idom[idom[rev[i]]];
        }
    }
};
