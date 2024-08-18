//区间加等差数列
void add(int l,int r,int a,int k){//a首项，k公差，d数组两次前缀和即为答案
    d2[l]+=a;
    d2[l+1]+=k-a;
    d2[r+1]-=(r-l+1)*k+a;
    d2[r+2]-=(l-r)*k-a;
}

//区间段合并

vector<pair<int,int>> merge_segment(vector<pair<int,int>> q){
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

//loading
