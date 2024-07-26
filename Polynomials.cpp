namespace Poly{
    using namespace std;
    #ifdef LOCAL
    #define debug(...) fprintf(stderr, ##__VA_ARGS__)
    #else
    #define debug(...) void(0)
    #endif
    typedef long long LL;
    template <unsigned M_> struct ModInt {
	    static constexpr unsigned M = M_;
	    unsigned x;
	    constexpr ModInt() : x(0U) {}
	    constexpr ModInt(unsigned x_) : x(x_ % M) {}
	    constexpr ModInt(unsigned long long x_) : x(x_ % M) {}
	    constexpr ModInt(int x_) : x(((x_ %= static_cast<int>(M)) < 0) ? (x_ + static_cast<int>(M)) : x_) {}
	    constexpr ModInt(long long x_) : x(((x_ %= static_cast<long long>(M)) < 0) ? (x_ + static_cast<long long>(M)) : x_) {}
	    ModInt(const std::string& str) {
	        x = 0;
	        size_t i = 0;
	        if (str.front() == '-') i += 1;
	        while (i < str.size()) {
	            assert(isdigit(str[i]));
	            x = (x * 10ull % M + str[i] - '0') % M;
	            i += 1;
	        }
	        if (str.front() == '-' && x) x = M - x;
	    }
	    ModInt &operator+=(const ModInt &a) { x = ((x += a.x) >= M) ? (x - M) : x; return *this; }
	    ModInt &operator-=(const ModInt &a) { x = ((x -= a.x) >= M) ? (x + M) : x; return *this; }
	    ModInt &operator*=(const ModInt &a) { x = (static_cast<unsigned long long>(x) * a.x) % M; return *this; }
	    ModInt &operator/=(const ModInt &a) { return (*this *= a.inv()); }
	    ModInt pow(long long e) const {
	        if (e < 0) return inv().pow(-e);
	        ModInt a = *this, b = 1U; for (; e; e >>= 1) { if (e & 1) b *= a; a *= a; } return b;
	    }
	    ModInt inv() const {
	        unsigned a = M, b = x; int y = 0, z = 1;
	        for (; b; ) { const unsigned q = a / b; const unsigned c = a - q * b; a = b; b = c; const int w = y - static_cast<int>(q) * z; y = z; z = w; }
	        assert(a == 1U); return ModInt(y);
	    }
	    ModInt operator+() const { return *this; }
	    ModInt operator-() const { ModInt a; a.x = x ? (M - x) : 0U; return a; }
	    ModInt operator+(const ModInt &a) const { return (ModInt(*this) += a); }
	    ModInt operator-(const ModInt &a) const { return (ModInt(*this) -= a); }
	    ModInt operator*(const ModInt &a) const { return (ModInt(*this) *= a); }
	    ModInt operator/(const ModInt &a) const { return (ModInt(*this) /= a); }
	    template <class T> friend ModInt operator+(T a, const ModInt &b) { return (ModInt(a) += b); }
	    template <class T> friend ModInt operator-(T a, const ModInt &b) { return (ModInt(a) -= b); }
	    template <class T> friend ModInt operator*(T a, const ModInt &b) { return (ModInt(a) *= b); }
	    template <class T> friend ModInt operator/(T a, const ModInt &b) { return (ModInt(a) /= b); }
	    explicit operator bool() const { return x; }
	    bool operator==(const ModInt &a) const { return (x == a.x); }
	    bool operator!=(const ModInt &a) const { return (x != a.x); }
	    bool operator<(const ModInt &a) const { return (x < a.x); }
	    bool operator>(const ModInt &a) const { return (x > a.x); }
	    bool operator<=(const ModInt &a) const { return (x <= a.x); }
	    bool operator>=(const ModInt &a) const { return (x >= a.x); }
	    friend int raw(const ModInt &a) { return a.x; }
	    friend ModInt qpow(ModInt a, long long b) {
	        ModInt r = 1;
	        for (; b; b >>= 1, a *= a)
	            if (b & 1) r *= a;
	        return r;
	    }
	    friend std::ostream &operator<<(std::ostream &os, const ModInt &a) { return os << a.x; }
	    friend std::istream &operator>>(std::istream &is, ModInt &a) {int v;is >> v;a = ModInt(v);return is;}
	};
	constexpr unsigned MO = 998244353;
	using Mint = ModInt<MO>;
    int glim(const int& x) { return 1 << (32 - __builtin_clz(x - 1)); }
    int bitctz(const int& x) { return __builtin_ctz(x); }
    struct poly : vector<Mint> {
        poly() {}
        explicit poly(int n) : vector<Mint>(n) {}
        poly(const vector<Mint>& vec) : vector<Mint>(vec) {}
        poly(initializer_list<Mint> il) : vector<Mint>(il) {}
        Mint operator()(const Mint& x) const;
        poly& cut(int lim);
        void ntt(int op);
    };
    void print(const poly& a) {
        for (size_t i = 0; i < a.size(); i++) debug("%d, ", raw(a[i]));
        debug("\n");
    }
    istream& operator>>(istream& is, poly& a) {
        for (auto& x : a) is >> x;
        return is;
    }
    ostream& operator<<(ostream& os, const poly& a) {
        bool flag = false;
        for (auto& x : a) {
            if (flag)
                os << " ";
            else
                flag = true;
            os << x;
        }
        return os;
    }
    Mint poly::operator()(const Mint& x) const {
        const auto& a = *this;
        Mint res = 0;
        for (int i = (int)a.size() - 1; i >= 0; i--) {
            res = res * x + a[i];
        }
        return res;
    }
    poly& poly::cut(int lim) {
        resize(lim);
        return *this;
    }
    void poly::ntt(int op) {
        static bool wns_flag = false;
        static vector<Mint> wns;
        if (!wns_flag) {
            wns_flag = true;
            for (int j = 1; j <= 23; j++) {
                wns.push_back(qpow(Mint(3), raw(Mint(-1)) >> j));
            }
        }
        vector<Mint>& a = *this;
        int n = a.size();
        for (int i = 1, r = 0; i < n; i++) {
            r ^= n - (1 << (bitctz(n) - bitctz(i) - 1));
            if (i < r) std::swap(a[i], a[r]);
        }
        vector<Mint> w(n);
        for (int k = 1, len = 2; len <= n; k <<= 1, len <<= 1) {
            Mint wn = wns[bitctz(k)];
            for (int i = raw(w[0] = 1); i < k; i++) w[i] = w[i - 1] * wn;
            for (int i = 0; i < n; i += len) {
                for (int j = 0; j < k; j++) {
                    Mint x = a[i + j], y = a[i + j + k] * w[j];
                    a[i + j] = x + y, a[i + j + k] = x - y;
                }
            }
        }
        if (op == -1) {
            Mint iz = Mint(1) / n;
            for (int i = 0; i < n; i++) a[i] *= iz;
            reverse(a.begin() + 1, a.end());
        }
    }
    poly concalc(int n, vector<poly> vec,
                            const function<Mint(vector<Mint>)>& func) {
        int lim = glim(n);
        int m = vec.size();
        for (auto& f : vec) f.resize(lim), f.ntt(1);
        vector<Mint> tmp(m);
        poly ret(lim);
        for (int i = 0; i < lim; i++) {
            for (int j = 0; j < m; j++) tmp[j] = vec[j][i];
            ret[i] = func(tmp);
        }
        ret.ntt(-1);
        return ret;
    }
    poly getInv(const poly& a, int lim) {
        poly b{1 / a[0]};
        for (int len = 2; len <= glim(lim); len <<= 1) {
            poly c = vector<Mint>(a.begin(), a.begin() + min(len, (int)a.size()));
            b = concalc(len << 1, {b, c}, [](vector<Mint> vec) {
                        return vec[0] * (2 - vec[0] * vec[1]);
                    }).cut(len);
        }
        return b.cut(lim);
    }
    poly operator+=(poly& a, const poly& b) {
        if (a.size() < b.size()) a.resize(b.size());
        for (size_t i = 0; i < b.size(); i++) a[i] += b[i];
        return a;
    }
    poly operator-=(poly& a, const poly& b) {
        if (a.size() < b.size()) a.resize(b.size());
        for (size_t i = 0; i < b.size(); i++) a[i] -= b[i];
        return a;
    }
    poly operator*=(poly& a, const Mint& k) {
        if (k == 1) return a;
        for (size_t i = 0; i < a.size(); i++) a[i] *= k;
        return a;
    }
    poly operator/=(poly& a, const Mint& k) { return a *= 1 / k; }
    poly operator<<=(poly& a, const int& k) {
        // mnltiple by x^k
        a.insert(a.begin(), k, 0);
        return a;
    }
    poly operator>>=(poly& a, const int& k) {
        // divide by x^k
        a.erase(a.begin(), a.begin() + min(k, (int)a.size()));
        return a;
    }
    poly operator*(const poly& a, const poly& b) {
        if (a.empty() || b.empty()) return {};
        int rlen = a.size() + b.size() - 1;
        int len = glim(rlen);
        if (1ull * a.size() * b.size() <= 1ull * len * bitctz(len)) {
            poly ret(rlen);
            for (size_t i = 0; i < a.size(); i++)
                for (size_t j = 0; j < b.size(); j++) ret[i + j] += a[i] * b[j];
            return ret;
        } else {
            return concalc(len, {a, b},
                                        [](vector<Mint> vec) { return vec[0] * vec[1]; })
                    .cut(rlen);
        }
    }
    poly operator/(poly a, poly b) {
        if (a.size() < b.size()) return {};
        int rlen = a.size() - b.size() + 1;
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());
        a = (a * getInv(b, rlen)).cut(rlen);
        reverse(a.begin(), a.end());
        return a;
    }
    poly operator-(poly a, const poly& b) { return a -= b; }
    poly operator%(const poly& a, const poly& b) {
        return (a - (a / b) * b).cut(b.size() - 1);
    }
    poly operator*=(poly& a, const poly& b) { return a = a * b; }
    poly operator/=(poly& a, const poly& b) { return a = a / b; }
    poly operator%=(poly& a, const poly& b) { return a = a % b; }
    poly operator+(poly a, const poly& b) { return a += b; }
    poly operator*(poly a, const Mint& k) { return a *= k; }
    poly operator*(const Mint& k, poly a) { return a *= k; }
    poly operator/(poly a, const Mint& k) { return a /= k; }
    poly operator<<(poly a, const int& k) { return a <<= k; }
    poly operator>>(poly a, const int& k) { return a >>= k; }
    poly getDev(poly a) {
        a >>= 1;
        for (int i = 1; i < (int)a.size(); i++) a[i] *= Mint(i) + 1;
        return a;
    }
    poly getInt(poly a) {
        a <<= 1;
        for (int i = 1; i < (int)a.size(); i++) a[i] /= Mint(i);
        return a;
    }
    poly getLn(const poly& a, int lim) {
        assert(a[0] == 1);
        return getInt(getDev(a) * getInv(a, lim)).cut(lim);
    }
    poly getExp(const poly& a, int lim) {
        assert(a[0] == 0);
        poly b{1};
        for (int len = 2; len <= glim(lim); len <<= 1) {
            poly c = vector<Mint>(a.begin(), a.begin() + min(len, (int)a.size()));
            b = concalc(len << 1, {b, getLn(b, len), c}, [](vector<Mint> vec) {
                        return vec[0] * (1 - vec[1] + vec[2]);
                    }).cut(len);
        }
        return b.cut(lim);
    }
    poly qpow(const poly& a, string k, int lim) {
        size_t i = 0;
        while (i < a.size() && a[i] == 0) i += 1;
        if (i == a.size() || (i > 0 && k.size() >= 9) ||
                1ull * i * raw(Mint(k)) >= 1ull * lim)
            return poly(lim);
        lim -= i * raw(Mint(k));
        return getExp(getLn(a / a[i] >> i, lim) * k, lim) *
                            qpow(a[i], raw(ModInt<Mint::M - 1>(k)))
                    << i * raw(Mint(k));
    }
    poly qpow(const poly& a, LL k, int lim) {
        size_t i = 0;
        while (i < a.size() && a[i] == 0) i += 1;
        if (i == a.size() || (i > 0 && k >= 1e9) ||
                1ull * i * k >= 1ull * lim)
            return poly(lim);
        lim -= i * k;
        return getExp(getLn(a / a[i] >> i, lim) * k, lim) *
                            qpow(a[i], raw(ModInt<Mint::M - 1>(k)))
                    << i * k;
    }
    Mint sqrt(const Mint& c) {
        static const auto check = [](Mint c) {
            return qpow(c, (Mint::M - 1) >> 1) == 1;
        };
        if (raw(c) <= 1) return 1;
        if (!check(c)) throw "No solution!";
        static mt19937 rng{random_device{}()};
        Mint a;a.x=rng();
        while (check(a * a - c)) {
            a.x = rng();
        }
        typedef pair<Mint, Mint> number;
        const auto mul = [=](number x, number y) {
            return make_pair(x.first * y.first + x.second * y.second * (a * a - c),
                                            x.first * y.second + x.second * y.first);
        };
        const auto qpow = [=](number a, int b) {
            number r = {1, 0};
            for (; b; b >>= 1, a = mul(a, a))
                if (b & 1) r = mul(r, a);
            return r;
        };
        Mint ret = qpow({a, 1}, (Mint::M + 1) >> 1).first;
        return min(raw(ret), raw(-ret));
    }
    poly getSqrt(const poly& a, int lim) {
        poly b{sqrt(a[0])};
        for (int len = 2; len <= glim(lim); len <<= 1) {
            poly c = vector<Mint>(a.begin(), a.begin() + min(len, (int)a.size()));
            b = (c * getInv(b * 2, len) + b / 2).cut(len);
        }
        return b.cut(lim);
    }
    template <class T>
    Mint divide_at(poly f, poly g, T n) {
        for (; n; n >>= 1) {
            poly r = g;
            for (size_t i = 1; i < r.size(); i += 2) r[i] *= -1;
            f *= r;
            g *= r;
            for (size_t i = n & 1; i < f.size(); i += 2) f[i >> 1] = f[i];
            f.resize((f.size() + 1) >> 1);
            for (size_t i = 0; i < g.size(); i += 2) g[i >> 1] = g[i];
            g.resize((g.size() + 1) >> 1);
        }
        return f.empty() ? 0 : f[0] / g[0];
    }
    template <class T>
    Mint linear_rec(poly a, poly f, T n) {
        // a[n] = sum_i f[i] * a[n - i]
        a.resize(f.size() - 1);
        f = poly{1} - f;
        poly g = a * f;
        g.resize(a.size());
        return divide_at(g, f, n);
    }
    poly BM(poly a) {
        poly ans, lst;
        int w = 0;
        Mint delta = 0;
        for (size_t i = 0; i < a.size(); i++) {
            Mint tmp = -a[i];
            for (size_t j = 0; j < ans.size(); j++) tmp += ans[j] * a[i - j - 1];
            if (tmp == 0) continue;
            if (ans.empty()) {
                w = i;
                delta = tmp;
                ans = vector<Mint>(i + 1, 0);
            } else {
                auto now = ans;
                Mint mul = -tmp / delta;
                if (ans.size() < lst.size() + i - w) ans.resize(lst.size() + i - w);
                ans[i - w - 1] -= mul;
                for (size_t j = 0; j < lst.size(); j++) ans[i - w + j] += lst[j] * mul;
                if (now.size() <= lst.size() + i - w) {
                    w = i;
                    lst = now;
                    delta = tmp;
                }
            }
        }
        return ans << 1;
    }
    poly lagrange(const vector<pair<Mint, Mint>>& a) {
        poly ans(a.size()), product{1};
        for (size_t i = 0; i < a.size(); i++) {
            product *= poly{-a[i].first, 1};
        }
        auto divide2 = [&](poly a, Mint b) {
            poly res(a.size() - 1);
            for (size_t i = (int)a.size() - 1; i >= 1; i--) {
                res[i - 1] = a[i];
                a[i - 1] -= a[i] * b;
            }
            return res;
        };
        for (size_t i = 0; i < a.size(); i++) {
            Mint denos = 1;
            for (size_t j = 0; j < a.size(); j++) {
                if (i != j) denos *= a[i].first - a[j].first;
            }
            poly numes = divide2(product, -a[i].first);
            ans += a[i].second / denos * numes;
        }
        return ans;
    }
}
