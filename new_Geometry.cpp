struct Point {
    double x, y;
    Point(double x0 = 0, double y0 = 0) : x(x0), y(y0) {}
    friend bool operator<(Point a, Point b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    }
    friend bool operator==(Point a, Point b) {
        return a.x == b.x and a.y==b.y;
    }
    friend Point operator+(const Point &a, const Point &b) {
        return Point(a.x + b.x, a.y + b.y);
    }
    friend Point operator-(const Point &a, const Point &b) {
        return Point(a.x - b.x, a.y - b.y);
    }
    friend Point operator*(const Point &a, double b) {
        return Point(a.x * b, a.y * b);
    }
    friend Point operator/(const Point &a, double b) {
        return Point(a.x / b, a.y / b);
    }
};

typedef Point Vector;

int dcmp(double x, double y) {
	if(fabs(x-y)<eps) return 0;
	if(x<y) return -1;
	return 1;
}

double Dot(const Point &a, const Point &b){
  	return a.x*b.x+a.y*b.y;
}

double Cross(const Point &a, const Point &b){
  	return a.x*b.y-a.y*b.x;
}

double Length(const Point &a){
  	return sqrt(Dot(a, a));
}

Point Normal(const Point &v,bool clockwise=1){
    double L=Length(v);
    if(!clockwise) return Point(-v.y/L,v.x/L);
    return Point(v.y/L,-v.x/L);
}

Point AcuteAngleFootPoint(const Point &P,const Point &A,const Point &B){
    Vector v1=P-A,v2=B-A;
    double len=Dot(v2,v1);
    double k=Length(v2);
    len/=k;
    return A+v2*(len/k);
}


bool PointOnLineLeft(const Point &P,const Point &A,const Point &B){
    Vector v=B-A,test=P-A;
    return dcmp(Cross(v,test),0)==1;
}

double PointDistToLine(const Point &P,const Point &A,const Point &B){
    Vector l=B-A,r=P-A;
    double v=Cross(l,r);
    return v/Length(l);
}

double PointDistToSegment(const Point &P,const Point &A,const Point &B){
	if(A==B) return Length(P-A);
	Vector v1=B-A,v2=P-A,v3=P-B;
	if(Dot(v1,v2)<0) return Length(v2);
	if(Dot(v1,v3)>0) return Length(v3);
	return fabs(PointDistToLine(P,A,B));
}

std::vector<int> ConvexHull(std::vector<Point> &p){
    int n=p.size();
    std::sort(p.begin(),p.end());
    std::vector<int> res;
    int m=0;
    for(int i=0;i<n;i++){
        while(m>=2 and dcmp(Cross(p[res[m-1]]-p[res[m-2]],p[i]-p[res[m-2]]),0)==-1){
            res.pop_back();
            m--;
        }
        res.push_back(i);m++;
    }
    int k=m;
    for(int i=n-2;i>=0;i--){
        while(m>k and dcmp(Cross(p[res[m-1]]-p[res[m-2]],p[i]-p[res[m-2]]),0)==-1){
            res.pop_back();
            m--;
        }
        res.push_back(i);m++;
    }
    res.pop_back();
    return res;
}

int MaxTriangleOnConvex(std::vector<Point> p,std::vector<int> cvx){
    if(cvx.size()<3) return 0;
    int m=cvx.size();
    int res=0;
    for(int i=0;i<m-1;i++){
        for(int j=i+1,k=(i+2)%m;j<m;j++){
            while(getDistToLine(p[cvx[k]],p[cvx[i]],p[cvx[j]])<getDistToLine(p[cvx[(k+1)%m]],p[cvx[i]],p[cvx[j]])) k=(k+1)%m;
            res=std::max(res,Cross(p[cvx[j]]-p[cvx[i]],p[cvx[k]]-p[cvx[i]]));
        }
    }
    return res;
}

int MinTriangleOnPoint(std::vector<Point> p){
    if(p.size()<3) return 0;
    int n=p.size();
    sort(p.begin(),p.end());
    std::vector<std::array<int,2>> K;
    for(int i=0;i<n;i++) {
        for(int j=i+1;j<n;j++) {
            if(i==j) continue;
            if(p[i]==p[j]) return 0;
            K.push_back({i,j});
        }
    }
    sort(K.begin(),K.end(),[&](auto l,auto r)->bool{
        if(p[l[0]].x==p[l[1]].x) return 1;
        if(p[r[0]].x==p[r[1]].x) return 0;
        int vl=(p[l[0]].y-p[l[1]].y)*(p[r[0]].x-p[r[1]].x);
        int vr=(p[l[0]].x-p[l[1]].x)*(p[r[0]].y-p[r[1]].y);
        return vl<vr;
    }); 
    std::vector<int> idx(n);
    std::vector<int> q(n);
    for(int i=0;i<n;i++) idx[i]=i,q[i]=i;

    int res=2e18;
    for(auto [i,j]:K){
        int l=idx[i],r=idx[j];
        if(l>0) res=std::min(res,std::abs(Cross(p[j]-p[i],p[q[l-1]]-p[i])));
        if(r+1<n) res=std::min(res,std::abs(Cross(p[j]-p[i],p[q[r+1]]-p[i])));
        if(p[i].x!=p[j].x){
            std::swap(q[l],q[r]);
            std::swap(idx[i],idx[j]);
        }
    }
    return res;
}

double MinRectangleCover(std::vector<Point> p,std::vector<int> cvx){
    double res=1e18;
    int n=cvx.size();
    assert((int)n>=3);
	std::vector<Point> ans(4);
    int J=2,L=2,R=1;
    for(int i=0;i<n;i++){
        int l=i,r=(i+1)%n;
        while(dcmp(PointDistToLine(p[cvx[J]],p[cvx[l]],p[cvx[r]]),PointDistToLine(p[cvx[(J+1)%n]],p[cvx[l]],p[cvx[r]]))==-1) J=(J+1)%n;
        Vector A=p[cvx[r]]-p[cvx[l]];
        while(1){
            Vector b=p[cvx[R]]-p[cvx[l]];
            Vector c=p[cvx[(R+1)%n]]-p[cvx[l]];
            if(dcmp(Dot(A,b),Dot(A,c))==-1){
                R=(R+1)%n;
            }else break;
        } 
        Vector B=p[cvx[l]]-p[cvx[r]];
        while(1){
            Vector b=p[cvx[L]]-p[cvx[r]];
            Vector c=p[cvx[(L+1)%n]]-p[cvx[r]];
            if(dcmp(Dot(B,b),Dot(B,c))!=1){
                L=(L+1)%n;
            }else break;
        }
        double h=PointDistToLine(p[cvx[J]],p[cvx[l]],p[cvx[r]]);
        double l1=Dot(p[cvx[R]]-p[cvx[l]],A)/Length(A),l2=Dot(p[cvx[L]]-p[cvx[r]],B)/Length(B);
        double len=l1+l2-Length(A);
        if(dcmp(len*h,res)==-1){
            res=len*h;
            ans[0]=AcuteAngleFootPoint(p[cvx[L]],p[cvx[l]],p[cvx[r]]);
            ans[1]=AcuteAngleFootPoint(p[cvx[R]],p[cvx[l]],p[cvx[r]]);
            ans[2]=ans[1]+Normal(ans[0]-ans[1])*h;
            ans[3]=ans[0]+Normal(ans[1]-ans[0],0)*h;
        }
    }

    return res;
}

double MinTwoConvexHullDist(std::vector<Point> p,std::vector<int> cvxp,std::vector<Point> q,std::vector<int> cvxq){
    double res=4e18;
    int n=cvxp.size(),m=cvxq.size();
    int idxp=0,idxq=0;
    for(int i=0;i<n;i++) if(dcmp(p[cvxp[i]].y,p[cvxp[idxp]].y)==-1) idxp=i;
    for(int i=0;i<m;i++) if(dcmp(q[cvxq[i]].y,q[cvxq[idxq]].y)==1) idxq=i;
 
    for(int k=0;k<n;k++){
        while(1){
            double l=Cross(q[cvxq[idxq]]-p[cvxp[idxp]],p[cvxp[(idxp+1)%n]]-p[cvxp[idxp]]);
            double r=Cross(q[cvxq[(idxq+1)%m]]-p[cvxp[idxp]],p[cvxp[(idxp+1)%n]]-p[cvxp[idxp]]);
            if(dcmp(l,r)!=-1) idxq=(idxq+1)%m;
            else break;
        }
         
        res=std::min(res,PointDistToSegment(q[cvxq[idxq]],p[cvxp[(idxp+1)%n]],p[cvxp[idxp]]));
        res=std::min(res,PointDistToSegment(q[cvxq[(idxq+1)%m]],p[cvxp[(idxp+1)%n]],p[cvxp[idxp]]));
        res=std::min(res,PointDistToSegment(p[cvxp[idxp]],q[cvxq[(idxq+1)%m]],q[cvxq[idxq]]));
        res=std::min(res,PointDistToSegment(p[cvxp[(idxp+1)%n]],q[cvxq[(idxq+1)%m]],q[cvxq[idxq]]));
 
        idxp=(idxp+1)%n;
    }
    return res;
}

struct Line{
    Point P;
    Vector v;
    double ang;
    Line(){}
    Line(Point P,Vector v):P(P),v(v){ 
        ang = atan2(v.y, v.x);
    }
    friend bool operator<(Line a,Line b){
        return a.ang < b.ang; 
    }
};
 
Point TwoLineIntersection(Line a,Line b){
    Vector u=a.P-b.P;
    double t=Cross(b.v,u)/Cross(a.v,b.v);
    return a.P+a.v*t;
}