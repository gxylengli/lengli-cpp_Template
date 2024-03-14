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