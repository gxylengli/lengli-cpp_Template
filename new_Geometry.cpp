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
};

typedef Point Vector;

double dcmp(double x, double y) {
  if (fabs(x - y) < eps) return 0;
  if (x < y) return -1;
  return 1;
}

double Dot(const Point &a, const Point &b){
  return a.x * b.x + a.y * b.y;
}

double Cross(const Point &a, const Point &b){
  return a.x * b.y - a.y * b.x;
}

double Length(const Point &a){
  return sqrt(Dot(a, a));
}

std::vector<int> ConvexHull(std::vector<Point> &p){
    int n=p.size();
    std::sort(p.begin(),p.end());
    std::vector<int> res;
    int m=0;
    for(int i=0;i<n;i++){
        while(m>=2 and Cross(p[res[m-1]]-p[res[m-2]],p[i]-p[res[m-2]])<=0){
            res.pop_back();
            m--;
        }
        res.push_back(i);m++;
    }
    int k=m;
    for(int i=n-2;i>=0;i--){
        while(m>k and Cross(p[res[m-1]]-p[res[m-2]],p[i]-p[res[m-2]])<=0){
            res.pop_back();
            m--;
        }
        res.push_back(i);m++;
    }
    res.pop_back();
    return res;
}