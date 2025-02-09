//区间加等差数列
void add(int l,int r,int a,int k){//a首项，k公差，d数组两次前缀和即为答案
    d2[l]+=a;
    d2[l+1]+=k-a;
    d2[r+1]-=(r-l+1)*k+a;
    d2[r+2]-=(l-r)*k-a;
}

//区间段合并

std::vector<pair<int,int>> merge_segment(vector<pair<int,int>> q){
    vector<pair<int,int>> p;
    sort(q.begin(),q.end());
    int l=-1,r=-1;
    for(pair<int,int> t:q){
        if(l==-1) l=t.first,r=t.second;
        else if(t.first>r) {
            p.push_back({l,r});
            t.first,r=t.second;
        }else r=max(r,t.second);
    }
    if(l!=-1) p.push_back({l,r});
    return p;
}

//字符串切分
std::vector<std::string> split(std::string a,std::string b){
	int n=a.size(),m=b.size();
	a="#"+a,b="#"+b;
	std::vector<int> q(n+1,0);
    std::vector<int> ne(m+1,0);
	for(int i=2,j=0;i<=m;i++){
		while(j and b[j+1]!=b[i]) j=ne[j];
		if(b[j+1]==b[i]) j++;
		ne[i]=j;
	}
	std::vector<std::string> res;
	std::string str="";
	for(int i=1,j=0;i<=n;i++){
		while(j and a[i]!=b[j+1]) j=ne[j];
		if(a[i]==b[j+1]) j++;
		str+=a[i];
		q[i]=j;
		if(j==m){
			j=0;
			str.erase(str.size()-m);
			if(str.size()) res.pb(str);
			str.clear();
		}
	}
	if(str.size()) res.pb(str);
	return res;
}


//动态合并树的直径
struct Diameter{
    int x,y;
    Diameter(){};
    Diameter(int x,int y):x(x),y(y){};
    friend Diameter operator + (Diameter a,Diameter b){
        int ma=0,rx=-1,ry=-1;
        for(auto l:{a.x,a.y,b.x,b.y}){
            for(auto r:{a.x,a.y,b.x,b.y}){ 
                int d=hld.dist(l,r);
                if(d>ma){
                    ma=d;
                    rx=l,ry=r;
                }
            }
        }
        return Diameter(rx,ry);
    }
};

//树的直径

int get_tree_diameter(int n,std::vector<std::vector<int>>& eg){
    std::vector<int> d(n+2,0);
    int c=0;
    auto dfs=[&](auto self,int u, int fa)->void{
        for(auto v:eg[u]) {
            if(v==fa) continue;
            d[v]=d[u]+1;
            if(d[v]>d[c]) c=v;
            self(self,v,u);
        }
    };
    dfs(dfs,1,0);
    d[c]=0;
    dfs(dfs,c,0);
    return d[c]; 
}

//unordered_map防止被卡

unordered_map<int, int> mp;
mp.reserve(1024);
mp.max_load_factor(0.25);

//矩阵旋转

std::vector<std::vector<int>> rotate(std::vector<std::vector<int>> &grid){
    int n=grid.size(),m=grid[0].size();
    std::vector<std::vector<int>> res(m,std::vector<int> (n));
    for(int i=0,jj=n-1;i<n;i++,jj--){
        for(int j=0,ii=0;j<m;j++,ii++){
            res[ii][jj]=grid[i][j];
        }
    }
    return res;
};

//手写哈希表

class Hash {
    private:
        int keys[N];
        int values[N];
    public:
        Hash() { memset(values, 0, sizeof(values)); }
        int& operator[](int n) {
            int idx = (n % N + N) % N, cnt = 1;
            while (keys[idx] != n and values[idx] != 0) {
                idx = (idx + cnt * cnt) % N;
                cnt += 1;
            }
            keys[idx] = n;
            return values[idx];
        }
}la;

//名次树，可以查排名的set

#include <bits/extc++.h>
using namespace __gnu_cxx;
using namespace __gnu_pbds;
using kth_tree = __gnu_pbds::tree<std::array<int,2>, null_type, std::less<>, rb_tree_tag, tree_order_statistics_node_update>;

//求逆序对
template<typename T>
long long inversion(std::vector<T> a,T max_val){
    long long res=0;
    if(max_val<=(int)1e6){
        BIT tr(max_val+2);
        reverse(all(a));
        for(auto x:a){
            res+=tr.sum(x-1);
            tr.add(x,1);
        }
    }else{
        std::vector<T> q;
        for(auto x:a) q.pb(x);
        q.erase(unique(all(q)),q.end());
        auto find=[&](int x){
            return lower_bound(all(q),x)-q.begin()+1;
        };
        BIT tr(q.size()+2);
        reverse(all(a));
        for(auto x:a){
            x=find(x);
            res+=tr.sum(x-1);
            tr.add(x,1);
        }
    }
    return res;
}


//unorderd_set/map,手写哈希

namespace std {
    template<>
    struct hash<std::array<int,2>> {
        size_t operator()(const std::array<int,2>& s) const {
            return hash<int>()(s[0]) + hash<int>()(s[1]);
        }
    };
}

//找拓扑序中的位置唯一的点
//https://www.luogu.com.cn/problem/P11073

int get_only_topsort_point(int n,std::vector<std::array<int,2>> edge){
    std::vector<std::vector<int>> eg(n+2);
    std::vector<int> rd(n+2,0);
    for(auto [l,r]:edge) eg[l].pb(r),rd[r]++;

    std::vector<int> seq;
    std::queue<int> q;
    for(int i=1;i<=n;i++) if(!rd[i]) q.push(i);
    while(q.size()){
        auto t=q.front();
        q.pop();
        seq.pb(t);
        for(auto x:eg[t]){
            rd[x]--;
            if(!rd[x]) q.push(x);
        }
    }

    std::vector<std::vector<int>> st(2,std::vector<int> (n+2,0));
    {
        rd.clear(),rd.resize(n+2,0);
        reverse(all(seq));
        int cnt=0;
        for(auto u:seq){
            for(int v:eg[u]){
                if(!rd[v]) cnt--;
                rd[v]++;
            }
            cnt++;
            if(cnt==1) st[0][u]=1;
        }
    }
    {
        eg.clear(),eg.resize(n+2);
        rd.clear(),rd.resize(n+2,0);
        for(auto [r,l]:edge) eg[l].pb(r);
        int cnt=0;
        reverse(all(seq));
        for(auto u:seq){
            for(int v:eg[u]){
                if(!rd[v]) cnt--;
                rd[v]++;
            }
            cnt++;
            if(cnt==1) st[1][u]=1;
        }
    }
    
    int res=0;
    for(int i=1;i<=n;i++){
        if(st[0][i] and st[1][i]) res++;
    }
    return res;
}

//找拓扑序中的位置唯一的点（时间戳做法）

int get_only_topsort_point(int n,std::vector<std::array<int,2>> edge){
    std::vector<std::vector<int>> eg(n+2);
    std::vector<int> rd(n+2,0);
    for(auto [l,r]:edge) {
        assert(l!=r);
        if(l>r) std::swap(l,r);
        eg[l].pb(r),rd[r]++;
    }

    std::vector<int> min_time(n+2,0);
    int tot=0;
    for(int i=1;i<=n;i++){
        for(int j:eg[i]){
            min_time[j]=max(min_time[j],min_time[i]+1);
            tot=max(tot,min_time[j]);
        }
    }
    std::vector<int> max_time(n+2,tot);
    for(int i=n;i>=1;i--){
        for(int j:eg[i]){
            max_time[i]=min(max_time[i],max_time[j]-1);
        }
    }
    std::vector<int> b(n+2);
    for(int i=1;i<=n;i++) b[min_time[i]]++,b[max_time[i]+1]--;
    for(int i=1;i<=n;i++) b[i]+=b[i-1];

    std::vector<int> res(n+2,1);
    for (int i=1;i<=n;i++) {
        bool feas=true;
        if(min_time[i]!=max_time[i]) feas=false;
        int t=min_time[i];
        if(b[t]>1) feas=false;
        res[i]=feas;
    }

    int ans=0;
    for(int i=1;i<=n;i++) if(res[i]) ans++;
    return ans;
}

//毫秒级随机数

std::mt19937 rd(std::chrono::system_clock::now().time_since_epoch().count());

//二维矩阵哈希

template<typename mint>
struct hash2D{
    int n=0,m=0;
    std::vector<mint> p1,p2;
    std::vector<std::vector<mint>> h;
    hash2D(){};
    hash2D(std::vector<std::vector<char>> &a,mint base1,mint base2){
        n=a.size(),m=a[0].size();
        h.resize(n+2,std::vector<mint> (m+2,0));
        int mx=std::max(n,m);
        p1.resize(mx+2,1),p2.resize(mx+2,1);
        for(int i=1;i<=mx;i++){
            p1[i]=p1[i-1]*base1;
            p2[i]=p2[i-1]*base2;
        }
        for(int i=1;i<=n;i++) {
            for(int j=1;j<=m;j++){
                h[i][j]=h[i][j-1]*base2;
                h[i][j]=h[i][j]+h[i-1][j]*base1;
                h[i][j]=h[i][j]-base2*base1*h[i-1][j-1];
                h[i][j]=h[i][j]+a[i-1][j-1];
            }
        } 
    };

    mint query(int x_1,int y_1,int x_2,int y_2){
        mint res=h[x_2][y_2];
        res=res-h[x_1-1][y_2]*p1[x_2-x_1+1];
        res=res-h[x_2][y_1-1]*p2[y_2-y_1+1];
        res=res+p1[x_2-x_1+1]*h[x_1-1][y_1-1]*p2[y_2-y_1+1];
        return res;
    }
};

//loading
