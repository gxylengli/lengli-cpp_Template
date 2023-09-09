#include <bits/stdc++.h>
#define fastio             \
  ios::sync_with_stdio(0); \
  cin.tie(0);              \
  cout.tie(0)
#define endl '\n'
//#define x first
//#define y second

using namespace std;
typedef pair<int, int> PII;
typedef long long LL;

const int N = 100010;
const double pi = acos(-1), eps = 1e-9;

//<----------------------------------------->点
struct Point {
  double x, y;
  Point(double x0 = 0, double y0 = 0) : x(x0), y(y0) {}
  friend bool operator<(Point a, Point b) {
    return a.x != b.x ? a.x < b.x : a.y < b.y;
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

double dcmp(double x, double y)  //比较两个数字大小
{
  if (fabs(x - y) < eps) return 0;
  if (x < y) return -1;
  return 1;
}
double dcmp(double x)  //符号判断
{
  if (fabs(x) < eps)
    return 0;
  else
    return x < 0 ? -1 : 1;
}
double Dot(const Point &a, const Point &b)  //点积
{
  return a.x * b.x + a.y * b.y;
}
double Cross(const Point &a, const Point &b)  //叉积
{
  return a.x * b.y - a.y * b.x;
}
double Length(const Point &a)  //求向量的长度
{
  return sqrt(Dot(a, a));
}
double Angle(const Point &a, const Point &b)  //求两个向量的夹角（余弦定理）
{
  return acos(Dot(a, b) / Length(a) / Length(b));
}
Point Rotate(const Point &a, double rad)  //逆时针旋转rad
{
  return Point(a.x * cos(rad) - a.y * sin(rad),
               a.x * sin(rad) + a.y * cos(rad));
}
Point Normal(Point &v)  //求单位长度的法向量
{
  double L = Length(v);
  return Point(-v.y / L, v.x / L);
}
Point getLineIntersection(const Point &P, const Point &v, const Point &Q,
                          const Point &w)  //求两个线段交点
{
  Vector u = P - Q;
  double t = Cross(w, u) / Cross(v, w);
  return P + v * t;
}
bool SegmentProperIntersection(
    const Point &a1, const Point &b1, const Point &a2,
    const Point &b2)  //两线段规范相交、即每条线段的端点分别在另一条一段的两侧
{
  double c1 = Cross(b1 - a1, a2 - a1), c2 = Cross(b1 - a1, b2 - a1);
  double c3 = Cross(b2 - a2, a1 - a2), c4 = Cross(b2 - a2, b1 - a2);
  return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
}
bool IsPointOnSegment(Point p, Point a1, Point a2)  //点在线段上
{
  return dcmp(Cross(p - a1, p - a2)) == 0 && dcmp(Dot(p - a1, p - a2)) < 0;
}
double getDistToLine(const Point &P, const Point &A,
                     const Point &B)  //点到直线的有向距离（距离加绝对值）
{
  Vector AB = B - A, AP = P - A;
  return Cross(AB, AP) / Length(AB);
}
int ConvexHull(Point *p, int n, Point *ch)  //构造逆时针凸包
{
  sort(p + 1, p + n + 1);  //先按照横坐标再按照纵坐标排序
  int m = 0;
  for (int i = 1; i <= n; i++) {
    while (m > 1 && Cross(ch[m] - ch[m - 1], p[i] - ch[m - 1]) <= 0) m--;
    ch[++m] = p[i];
  }
  int k = m;
  for (int i = n - 1; i; i--) {
    while (m > k && Cross(ch[m] - ch[m - 1], p[i] - ch[m - 1]) <= 0) m--;
    ch[++m] = p[i];
  }
  if (n > 1) m--;
  return m;
}
double PolygonArea(Point *p, int n)  //求逆时针构成的多边形（可不凸）面积
{
  double ret = 0;
  for (int i = 2; i < n; i++)  //第一个点是p[1],所以这样循环
    ret += Cross(p[i] - p[1], p[i + 1] - p[1]);
  return ret / 2;
}
bool isPointInPolygon(Point p, Point *poly, int n)  //点在凸多边形内的判定
{
  int wn = 0;
  poly[n + 1] = poly[1];
  for (int i = 1; i <= n; i++) {
    if (IsPointOnSegment(p, poly[i], poly[i + 1])) return -1;
    int k = dcmp(Cross(poly[i + 1] - poly[i], p - poly[i]));
    int d1 = dcmp(poly[i].y - p.y);
    int d2 = dcmp(poly[i + 1].y - p.y);
    if (k > 0 && d1 <= 0 && d2 > 0) wn++;
    if (k < 0 && d2 <= 0 && d1 > 0) wn--;
  }
  if (wn != 0) return 1;
  return 0;
}
void update(int a, int b) {}
int findDiameter(Point *p, int n)  //旋转卡壳求直径（Diatmeter：直径）
{
  int mx = 0, j = 2;
  p[n + 1] = p[1];
  for (int i = 1; i <= n; i++) {
    while (dcmp(Cross(p[i + 1] - p[i], p[j + 1] - p[j])) > 0) {
      j++;
      j = (j - 1) % n + 1;
    }
    update(i, j);
    update(i + 1, j);
    update(i, j + 1);
    update(i + 1, j + 1);
  }
  return mx;
}
//<----------------------------------------->线
struct Line {
  Point P;   //起点
  Vector v;  //从起点到终点的向量
  double ang;
  Line() {}
  Line(Point P, Vector v) : P(P), v(v) { ang = atan2(v.y, v.x); }
  friend bool operator<(Line a, Line b) { return a.ang < b.ang; }
};

Point GetIntersection(Line a, Line b)  //两条直线交点
{
  Vector u = a.P - b.P;
  double t = Cross(b.v, u) / Cross(a.v, b.v);
  return a.P + a.v * t;
}
bool OnLeft(Line L, Point p) { return Cross(L.v, p - L.P) >= 0; }
int HalfplaneIntersection(Line *L, int n, Point *poly)  //半平面交
{
  sort(L + 1, L + n + 1);
  int first, last;
  Point *p = new Point[n + 10];
  Line *q = new Line[n + 10];
  q[first = last = 0] = L[1];
  for (int i = 2; i <= n; i++) {
    while (first < last && !OnLeft(L[i], p[last - 1])) last--;
    while (first < last && !OnLeft(L[i], p[first])) first++;
    q[++last] = L[i];
    if (fabs(Cross(q[last].v, q[last - 1].v)) < eps) {
      last--;
      if (OnLeft(q[last], L[i].P)) q[last] = L[i];
    }
    if (first < last) p[last - 1] = GetIntersection(q[last - 1], q[last]);
  }
  while (first < last && !OnLeft(q[first], p[last - 1])) last--;
  if (last - first <= 1) return 0;
  p[last] = GetIntersection(q[last], q[first]);
  int m = 0;
  for (int i = first; i <= last; i++) poly[++m] = p[i];
  return m;
}
//<----------------------------------------->圆
struct Circle {
  Point P;
  double r;
  Circle(){};
  Circle(Point a, double r) : P(a), r(r){};
};

double getPointdist(const Point &a, const Point &b)  //两点距离
{
  double dx = b.x - a.x;
  double dy = b.y - a.y;
  double len = sqrt(dx * dx + dy * dy);
  return len;
}

bool PointinCircle(const Point &a, const Circle &b)  //点是否在圆内
{
  double len = getPointdist(a, b.P);
  if (len < b.r)
    return 1;
  else
    return 0;
}

Circle getCircle(const Point &a, const Point &b, const Point &c)  //三点确定圆
{
  Point aa = b - a, bb = c - a;
  aa = Rotate(aa, -pi / 2);
  bb = Rotate(bb, -pi / 2);
  Line p = Line((b + a) / 2, aa);
  Line q = Line((c + a) / 2, bb);
  Point ans = GetIntersection(p, q);
  double r = getPointdist(ans, a);
  return Circle(ans, r);
}

Circle minCirclecover(Point *q, int n)  //最小圆覆盖
{
  random_shuffle(q + 1, q + 1 + n);
  Circle c = Circle(q[1], 0);
  for (int i = 2; i <= n; i++) {
    if (!PointinCircle(q[i], c)) {
      c = Circle(q[i], 0);
      for (int j = 1; j < i; j++) {
        if (!PointinCircle(q[j], c)) {
          c = Circle((q[i] + q[j]) / 2, getPointdist(q[i], q[j]) / 2);

          for (int k = 1; k < j; k++) {
            if (!PointinCircle(q[k], c)) {
              c = getCircle(q[i], q[j], q[k]);
            }
          }
        }
      }
    }
  }
  return c;
}
double AreaOfOverlap(Point c1, double r1, Point c2, double r2)  //两圆相交面积
{
  double d = Length(c1 - c2);
  if (r1 + r2 < d + eps) return 0.0;
  if (d < fabs(r1 - r2) + eps) {
    double r = min(r1, r2);
    return pi * r * r;
  }
  double x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
  double p = (r1 + r2 + d) / 2.0;
  double t1 = acos(x / r1);
  double t2 = acos((d - x) / r2);
  double s1 = r1 * r1 * t1;
  double s2 = r2 * r2 * t2;
  double s3 = 2 * sqrt(p * (p - r1) * (p - r2) * (p - d));
  return s1 + s2 - s3;
}
signed main() {
  fastio;

  return 0;
}