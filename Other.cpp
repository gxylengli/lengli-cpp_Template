//区间加等差数列
void add(int l,int r,int a,int k){//a首项，k公差，d数组两次前缀和即为答案
    d2[l]+=a;
    d2[l+1]+=k-a;
    d2[r+1]-=(r-l+1)*k+a;
    d2[r+2]-=(l-r)*k-a;
}