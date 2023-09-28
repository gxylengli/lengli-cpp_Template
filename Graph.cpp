//最短路dijkstra

vector<int> dijkstra(int S){
	vector<int> dist(n+1,-1);
	vector<bool> st(n+1,0);
	priority_queue<PII,vector<PII>,greater<PII>> q;
	q.push({0,S});dist[S]=0;
	while(q.size()){
		auto t=q.top();
		q.pop();
		if(st[t.second]) continue;
		st[t.second]=1;
		for(auto x:eg[t.second]){
			if(dist[x]==-1 or dist[x]>t.first+1){
				dist[x]=t.first+1;
				q.push({dist[x],x});
			}
		}
	}
	return dist;
}

//次短路dijkstra

void dijkstra2(int u) {
  for (int i = 0; i < N; i++) dis1[i] = dis2[i] = 1e18;
  priority_queue<PII, vector<PII>, greater<PII>> q;
  q.push({0, u});
  dis1[u] = 0;
  while (!q.empty()) {
    auto t = q.top();
    q.pop();
    if (dis2[t.second] < t.first) continue;
    for (auto [w, j] : eg[t.second]) {
      int nw = t.first + w;
      if (nw < dis1[j]) {
        swap(dis1[j], nw);
        q.push({dis1[j], j});
      }
      if (nw < dis2[j]) {
        swap(dis2[j], nw);
        q.push({dis2[j], j});
      }
    }
  }
}

//匈牙利

#include <algorithm>
#include <cstring>
#include <iostream>

using namespace std;

const int N = 510, M = 100010;

int n1, n2, m;
int h[N], e[M], ne[M], idx;
int match[N];
bool st[N];

bool find(int x) {
  for (int i = h[x]; i != -1; i = ne[i]) {
    int j = e[i];
    if (!st[j]) {
      st[j] = true;
      if (match[j] == 0 || find(match[j])) {
        match[j] = x;
        return true;
      }
    }
  }
  return false;
}
void add(int a, int b) { e[idx] = b, ne[idx] = h[a], h[a] = idx++; }

int main() {
  cin >> n1 >> n2 >> m;
  memset(h, -1, sizeof h);
  while (m--) {
    int a, b;
    cin >> a >> b;
    add(a, b);
  }
  int res = 0;
  for (int i = 1; i <= n1; i++) {
    memset(st, false, sizeof st);
    if (find(i)) res++;
  }

  cout << res;

  return 0;
}

//最小生成树（kruskal）

#include <algorithm>
#include <iostream>

using namespace std;

const int N = 200010;

int n, m;
int p[N];

struct Edge {
  int a, b, w;

  bool operator<(const Edge &x) const { return w < x.w; }
} edges[N];

int find(int x) {
  if (p[x] != x) p[x] = find(p[x]);
  return p[x];
}

int main() {
  cin >> n >> m;
  for (int i = 0; i < m; i++) {
    int a, b, w;
    cin >> a >> b >> w;
    edges[i] = {a, b, w};
  }

  sort(edges, edges + m);

  for (int i = 1; i <= n; i++) p[i] = i;

  int res = 0, cnt = 0;

  for (int i = 0; i < m; i++) {
    int a = edges[i].a, b = edges[i].b, w = edges[i].w;

    a = find(a), b = find(b);

    if (a != b) {
      p[a] = b;
      cnt++;
      res += w;
    }
  }
  if (cnt < n - 1)
    cout << "impossible";
  else
    cout << res;
}

//倍增LCA（最近公共祖先）

#include <bits/stdc++.h>
#define int long long

using namespace std;
const int N = 4e4 + 10, M = 2 * N;

int n, m, root;
int h[N], e[M], ne[M], idx;
int dep[N];
int fa[N][16];

void add(int a, int b) { e[idx] = b, ne[idx] = h[a], h[a] = idx++; }

void dfs(int u, int f) {
  for (int i = 1; i < 16; i++) {
    fa[u][i] = fa[fa[u][i - 1]][i - 1];
  }
  for (int i = h[u]; ~i; i = ne[i]) {
    int j = e[i];
    if (j == f) continue;
    fa[j][0] = u;
    dep[j] = dep[u] + 1;
    dfs(j, u);
  }
}

int lca(int a, int b) {
  if (dep[a] < dep[b]) swap(a, b);
  for (int i = 15; i >= 0; i--) {
    if (dep[a] - (1 << i) >= dep[b]) a = fa[a][i];
  }
  if (a == b) return a;
  for (int i = 15; i >= 0; i--) {
    if (fa[a][i] != fa[b][i]) {
      a = fa[a][i];
      b = fa[b][i];
    }
  }
  return fa[a][0];
}

signed main() {
  ios::sync_with_stdio(0);
  cin >> n;
  memset(h, -1, sizeof h);
  for (int i = 0; i < n; i++) {
    int a, b;
    cin >> a >> b;
    if (!~b)
      root = a;
    else
      add(a, b), add(b, a);
  }
  dep[root] = 0;
  dfs(root, -1);

  cin >> m;
  while (m--) {
    int a, b;
    cin >> a >> b;
    int res = lca(a, b);

    if (res == a)
      cout << 1 << endl;
    else if (res == b)
      cout << 2 << endl;
    else
      cout << 0 << endl;
  }
}

//树链剖分，修改路径，修改子树
#include <bits/stdc++.h>
#define fastio             \
  ios::sync_with_stdio(0); \
  cin.tie(0);              \
  cout.tie(0)
#pragma GCC optimize(2)
#define endl '\n'
//#define x first
//#define y second

using namespace std;
typedef pair<int, int> PII;
typedef long long LL;

const int N = 100010;

int n, m;
int a[N];
vector<int> eg[N];
int dep[N], son[N], fa[N], sz[N];
int dfn[N], cnt, w[N], id[N], top[N];

void dfs1(int u, int f, int d) {
  dep[u] = d;
  sz[u] = 1;
  fa[u] = f;
  for (auto x : eg[u]) {
    if (x == f) continue;
    dfs1(x, u, d + 1);
    sz[u] += sz[x];
    if (sz[x] > sz[son[u]]) son[u] = x;
  }
}
void dfs2(int u, int f) {
  dfn[u] = ++cnt, w[cnt] = a[u];
  top[u] = f;
  if (!son[u]) return;
  dfs2(son[u], f);
  for (auto x : eg[u]) {
    if (x == fa[u] or x == son[u]) continue;
    dfs2(x, x);
  }
}

struct Node {
  int l, r;
  LL sum, add;
} tr[N * 4];
void pushup(Node &u, Node &l, Node &r) { u.sum = l.sum + r.sum; }
void pushup(int u) { pushup(tr[u], tr[u << 1], tr[u << 1 | 1]); }
void build(int u, int l, int r) {
  tr[u] = {l, r, 0, 0};
  if (l == r)
    tr[u] = {l, r, w[l], 0};
  else {
    int mid = (l + r) / 2;
    build(u << 1, l, mid);
    build(u << 1 | 1, mid + 1, r);
    pushup(u);
  }
}
void pushdown(int u) {
  if (tr[u].add) {
    tr[u << 1].sum += (tr[u << 1].r - tr[u << 1].l + 1) * tr[u].add,
        tr[u << 1].add += tr[u].add;
    tr[u << 1 | 1].sum += (tr[u << 1 | 1].r - tr[u << 1 | 1].l + 1) * tr[u].add,
        tr[u << 1 | 1].add += tr[u].add;
    tr[u].add = 0;
  }
}
void modify(int u, int l, int r, int d) {
  if (tr[u].l >= l and tr[u].r <= r) {
    tr[u].sum += (tr[u].r - tr[u].l + 1) * d;
    tr[u].add += d;
    return;
  }
  pushdown(u);
  int mid = (tr[u].l + tr[u].r) / 2;
  if (l <= mid) modify(u << 1, l, r, d);
  if (r > mid) modify(u << 1 | 1, l, r, d);
  pushup(u);
}
Node query(int u, int l, int r) {
  if (tr[u].l >= l and tr[u].r <= r) return tr[u];
  pushdown(u);
  int mid = (tr[u].r + tr[u].l) / 2;
  if (l > mid)
    return query(u << 1 | 1, l, r);
  else if (r <= mid)
    return query(u << 1, l, r);
  else {
    Node ll = query(u << 1, l, r), rr = query(u << 1 | 1, l, r);
    Node res;
    pushup(res, ll, rr);
    return res;
  }
}

void modify_tree(int u, int k) { modify(1, dfn[u], dfn[u] + sz[u] - 1, k); }
LL query_tree(int u) { return query(1, dfn[u], dfn[u] + sz[u] - 1).sum; }
void modify_path(int u, int v, int k) {
  while (top[u] != top[v]) {
    if (dep[top[u]] < dep[top[v]]) swap(u, v);
    modify(1, dfn[top[u]], dfn[u], k);
    u = fa[top[u]];
  }
  if (dep[u] < dep[v]) swap(u, v);
  modify(1, dfn[v], dfn[u], k);
}
LL query_path(int u, int v) {
  LL res = 0;
  while (top[u] != top[v]) {
    if (dep[top[u]] < dep[top[v]]) swap(u, v);
    res += query(1, dfn[top[u]], dfn[u]).sum;
    u = fa[top[u]];
  }
  if (dep[u] < dep[v]) swap(u, v);
  res += query(1, dfn[v], dfn[u]).sum;
  return res;
}

void solve() {
  cin >> n;
  for (int i = 1; i <= n; i++) cin >> a[i];
  for (int i = 1; i < n; i++) {
    int l, r;
    cin >> l >> r;
    eg[l].emplace_back(r);
    eg[r].emplace_back(l);
  }
  dfs1(1, -1, 1);
  dfs2(1, 1);
  build(1, 1, n);

  cin >> m;
  while (m--) {
    int x;
    cin >> x;
    if (x == 1) {
      int u, v, k;
      cin >> u >> v >> k;
      modify_path(u, v, k);
    } else if (x == 2) {
      int u, k;
      cin >> u >> k;
      modify_tree(u, k);
    } else if (x == 3) {
      int u, v;
      cin >> u >> v;
      cout << query_path(u, v) << endl;
    } else {
      int u;
      cin >> u;
      cout << query_tree(u) << endl;
    }
  }
}

signed main() {
  fastio;

  int T;
  T = 1;
  while (T--) solve();

  return 0;
}

//最大流:节点编号(1-n)

struct MF {
  int h[N], e[M], ne[M], w[M], idx = 0;

  int n, S, T;
  LL maxflow = 0;
  int dep[N], cur[N];

  void init() {
    memset(h, -1, sizeof h);
    idx = 0;
  }

  void add(int a, int b, int c) {
    e[idx] = b, ne[idx] = h[a], w[idx] = c, h[a] = idx++;
    e[idx] = a, ne[idx] = h[b], w[idx] = 0, h[b] = idx++;
  }

  bool bfs() {
    queue<int> q;
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
        int t = dfs(j, min(limit - flow, w[i]));
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
      // memcpy(cur, h, sizeof h);
      maxflow += dfs(S, INF);
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
	
	inline void add(int x,int y,int z,int c,bool o=1){
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
		queue<int> q;
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
			k=dfs(y,min(flow,z));
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

// Tarjan求点双连通分量V-DCC

struct Tarjan{
    vector<int> eg[N];
    int dfn[N],low[N],timespace;
    int stk[N],top;
    bool cut[N];
    int dcc_cnt,n;
    vector<int> dcc[N];
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
                low[u]=min(low[u],low[j]);
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
            else low[u]=min(low[u],dfn[j]);
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
    vector<int> eg[N];
    int dfn[N],low[N],timespace;
    int stk[N],top;
    bool st[N];
    int color[N];
    int scc_cnt,n;
    vector<int> scc[N];
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
				low[u]=min(low[u],low[v]);
			}else if(st[v]) low[u]=min(low[u],low[v]);
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

