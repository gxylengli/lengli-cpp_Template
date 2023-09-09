#include <bits/stdc++.h>
#define fastio             \
  ios::sync_with_stdio(0); \
  cin.tie(0);
#define endl '\n'

using namespace std;

bool st[10];
vector<pair<int, int>> a, b, all;

// 5张牌得到rank
vector<int> getRank(vector<pair<int, int>> vec) {
  bool isFlush = true;  // first  值  second 花色
  // 判断是不是同花
  for (int i = 1; i < 5; i++) {
    if (vec[i].second != vec[i - 1].second) {
      isFlush = false;
      break;
    }
  }
  // 数字从大到小排列
  int a[5];
  for (int i = 0; i < 5; i++) a[i] = vec[i].first;
  sort(a, a + 5);
  reverse(a, a + 5);

  vector<int> ret;

  // 10 - 皇家同花顺
  if (isFlush && a[0] == 14 && a[1] == 13 && a[2] == 12 && a[3] == 11 &&
      a[4] == 10) {
    ret.push_back(10);  // rank 具体的牌
    for (int i = 0; i < 5; i++) ret.push_back(a[i]);
    return ret;
  }
  // 判断是不是顺子
  bool isStraight = false;
  // 特判 A 2 3 4 5
  if (a[0] == 14 && a[1] == 5 && a[2] == 4 && a[3] == 3 && a[4] == 2) {
    isStraight = true;
  }
  if (a[0] == a[1] + 1 && a[1] == a[2] + 1 && a[2] == a[3] + 1 &&
      a[3] == a[4] + 1) {
    isStraight = true;
  }

  // 9 - 同花顺
  if (isStraight && isFlush) {
    ret.push_back(9);
    if (a[0] == 14 && a[1] == 5 && a[2] == 4 && a[3] == 3 && a[4] == 2) {
      ret.push_back(5);
      ret.push_back(4);
      ret.push_back(3);
      ret.push_back(2);
      ret.push_back(1);
    } else {
      for (int i = 0; i < 5; i++) ret.push_back(a[i]);
    }
    return ret;
  }

  // 8 - 四条
  if (a[0] == a[3]) {
    ret.push_back(8);  // rank   0 1 2 3
    for (int i = 0; i < 4; i++) ret.push_back(a[i]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[1] == a[4]) {
    ret.push_back(8);
    for (int i = 1; i < 5; i++)  //先放对子
      ret.push_back(a[i]);
    ret.push_back(a[0]);
    return ret;
  }

  // 7 - 葫芦（三带二）
  if (a[0] == a[2] && a[3] == a[4]) {
    ret.push_back(7);
    for (int i = 0; i < 5; i++) ret.push_back(a[i]);
    return ret;
  }
  if (a[0] == a[1] && a[2] == a[4]) {
    ret.push_back(7);
    for (int i = 2; i < 5; i++) ret.push_back(a[i]);
    for (int i = 0; i < 2; i++) ret.push_back(a[i]);
    return ret;
  }

  // 6 - 同花
  if (isFlush) {
    ret.push_back(6);
    for (int i = 0; i < 5; i++) ret.push_back(a[i]);
    return ret;
  }

  // 5 - 顺子
  if (isStraight) {
    ret.push_back(5);
    if (a[0] == 14 && a[1] == 5 && a[2] == 4 && a[3] == 3 && a[4] == 2) {
      ret.push_back(5);
      ret.push_back(4);
      ret.push_back(3);
      ret.push_back(2);
      ret.push_back(1);
    } else {
      for (int i = 0; i < 5; i++) ret.push_back(a[i]);
    }
    return ret;
  }

  // 4 - 三条
  if (a[0] == a[2]) {
    ret.push_back(4);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[1] == a[3]) {
    ret.push_back(4);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[0]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[2] == a[4]) {
    ret.push_back(4);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    return ret;
  }

  // 3 - 两对
  if (a[0] == a[1] && a[2] == a[3]) {
    ret.push_back(3);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[1] == a[2] && a[3] == a[4]) {
    ret.push_back(3);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    ret.push_back(a[0]);
    return ret;
  }
  if (a[0] == a[1] && a[3] == a[4]) {
    ret.push_back(3);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    ret.push_back(a[2]);
    return ret;
  }

  // 2 - 一对
  if (a[0] == a[1]) {
    ret.push_back(2);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[1] == a[2]) {
    ret.push_back(2);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    ret.push_back(a[0]);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[2] == a[3]) {
    ret.push_back(2);
    ret.push_back(a[2]);
    ret.push_back(a[3]);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[4]);
    return ret;
  }
  if (a[3] == a[4]) {
    ret.push_back(2);
    ret.push_back(a[3]);
    ret.push_back(a[4]);
    ret.push_back(a[0]);
    ret.push_back(a[1]);
    ret.push_back(a[2]);
    return ret;
  }

  // 1 - 高牌
  ret.push_back(1);
  for (int i = 0; i < 5; i++) ret.push_back(a[i]);
  return ret;
}

// 比较两个vector大小(长度都为6)，字典序
int cmp(vector<int>& a, vector<int>& b) {
  for (int i = 0; i < 6; i++) {
    if (a[i] > b[i])
      return 1;
    else if (a[i] < b[i])
      return -1;
  }
  return 0;
}

pair<int, int> getCard(string str) {
  pair<int, int> ret;
  if (str[0] == 'A')
    ret.first = 14;
  else if (str[0] == 'K')
    ret.first = 13;
  else if (str[0] == 'Q')
    ret.first = 12;
  else if (str[0] == 'J')
    ret.first = 11;
  else if (str[0] == 'T')
    ret.first = 10;
  else
    ret.first = str[0] - '0';

  if (str[1] == 'S') ret.second = 1;
  if (str[1] == 'H') ret.second = 2;
  if (str[1] == 'C') ret.second = 3;
  if (str[1] == 'D') ret.second = 4;

  return ret;
}

int dfs(int step) {
  if (step == 7) {
    vector<int> aa = getRank(a);
    vector<int> bb = getRank(b);
    return cmp(aa, bb);
  }
  int Draw = 0;
  if (step & 1) {
    for (int i = 0; i < 6; i++) {
      if (st[i]) continue;
      st[i] = 1;
      a.push_back(all[i]);

      int ret = dfs(step + 1);

      st[i] = 0;
      a.pop_back();

      if (ret == 1)
        return 1;
      else if (ret == 0)
        Draw++;
    }
    if (Draw) return 0;
    return -1;
  } else {
    for (int i = 0; i < 6; i++) {
      if (st[i]) continue;
      st[i] = 1;
      b.push_back(all[i]);
      int ret = dfs(step + 1);

      st[i] = 0;
      b.pop_back();

      if (ret == -1)
        return -1;
      else if (ret == 0)
        Draw++;
    }
    if (Draw) return 0;
    return 1;
  }
}

void solve() {
  a = vector<pair<int, int>>();
  b = vector<pair<int, int>>();
  all = vector<pair<int, int>>();
  memset(st, 0, sizeof st);

  string s;
  cin >> s;
  a.push_back(getCard(s));
  cin >> s;
  a.push_back(getCard(s));
  cin >> s;
  b.push_back(getCard(s));
  cin >> s;
  b.push_back(getCard(s));

  for (int i = 0; i < 6; i++) {
    cin >> s;
    all.push_back(getCard(s));
  }

  int res = dfs(1);
  if (res == 1)
    cout << "Alice" << endl;
  else if (res == -1)
    cout << "Bob" << endl;
  else
    cout << "Draw" << endl;
}

signed main() {
  int T;
  cin >> T;
  while (T--) solve();

  return 0;
}

/*
9
JC 4H
TS 5D
JS JH JD 4S 4C 4D
JC 4H
TS 5D
TH TC TD 5S 5H 5C
JC 4H
TS 5D
4S JS 5S TH TC TD
7C 3C
7H TH
3S 3H 3D 2C 4H 5S
7C 3C
7H TH
2H 4H 5H 6H 8H 9H
7C 3C
7H TH
TS 3S 2S 2H 4C 4D
2D KH
4D JC
2S 2H 2C KS KC KD
2D KH
4D JC
4S 4H 4C JS JH JD
2D KH
4D JC
2S KS 4S JS JH JD
*/
//----
/*
Alice
Bob
Draw
Alice
Bob
Draw
Alice
Bob
Draw
*/
