template<typename T> struct matrix {
	int n, m;
	std::vector<T> data;
	matrix () : matrix(0, 0) {}
	matrix (int n) : matrix(n, n) {}
	matrix (int n, int m) : n(n), m(m), data(n * m) {}
	matrix (const std::vector<std::vector<T> > &a) : n(a.size()), m(a[0].size()) {
		data.resize(n*m);
		for (int i = 0; i < n; i++) {
			copy(a[i].begin(), a[i].end(), data.begin()+i*m);
		}
	}
	auto operator [] (int i) { return data.begin() + i * m; }
	auto operator [] (int i) const { return data.cbegin() + i * m; }
	static matrix mul_ident(int n) {
		matrix res(n, n);
		for (int i = 0; i < n; i++) res[i][i] = 1;
		return res;
	}
	matrix operator + (const matrix &rhs) {
		assert(n == rhs.n && m == rhs.m);
		matrix res(n,m);
		for (int i = 0; i < n ; i++) {
			for (int j = 0; j < m; j++)
				res[i][j] = (*this)[i][j] + rhs[i][j];
		}
		return res;
	}
	matrix operator - (const matrix &rhs) {
		assert(n == rhs.n && m == rhs.m);
		matrix res(n,m);
		for (int i = 0; i < n ; i++) {
			for (int j = 0; j < m; j++)
				res[i][j] = (*this)[i][j] - rhs[i][j];
		}
		return res;
	}
	matrix operator * (const matrix &rhs) {
		assert(m == rhs.n);
		matrix res(n, rhs.m);
		for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) for (int k = 0; k < rhs.m; k++)
			res[i][k] += (*this)[i][j] * rhs[j][k];
		return res;
	}
	template<typename int_t> matrix operator ^ (int_t x) {
		assert(n == m);
		matrix res = mul_ident(n);
        matrix a = *this;
		while (x) {
			if (x & 1) res = a * res;
			a = a * a;
			x>>=1;
		}
		return res;
	}
	matrix &operator += (const matrix &rhs) { return *this = *this + rhs;}
	matrix &operator -= (const matrix &rhs) { return *this = *this - rhs;}
	matrix &operator *= (const matrix &rhs) { return *this = *this * rhs;}
	template<typename int_t> matrix &operator ^= (int_t x) { return *this = *this ^ x;}
	bool operator == (const matrix &rhs) const { return m == rhs.m && data == rhs.data; }
	bool operator != (const matrix &rhs) const { return m != rhs.m || data != rhs.data; }
	std::pair<bool, matrix> inv() {
		assert(n == m);
		matrix a = *this;
		matrix r = mul_ident(n);
		for (int i = 0; i < n; i++) {
			int id = -1;
			for (int j = i; j < n; j++) if (a[j][i] != T(0)) { id = j; break; }
			if (id == -1) return {false, matrix()};
			for (int j = i; j < n; j++) std::swap(a[i][j], a[id][j]);
			for (int j = 0; j < n; j++) std::swap(r[i][j], r[id][j]);
			auto t = T(1) / a[i][i];
			for (int j = i; j < n; j++) a[i][j] *= t;
			for (int j = 0; j < n; j++) r[i][j] *= t;
			for (int j = 0; j < n; j++) if (i != j) {
				auto s = a[j][i];
				for (int k = i; k < n; k++) a[j][k] -= a[i][k] * s;
				for (int k = 0; k < n; k++) r[j][k] -= r[i][k] * s;
			}
		}
		return {true, r};
	}
	std::pair<bool, matrix> inv2() {
		assert(n == m);
		matrix a = *this;
		std::vector<std::pair<int, int> > swaps;
		for (int i = 0; i < n; i++) {
			int id = -1;
			for (int j = i; j < n; j++) if (a[j][i] != T(0)) { id = j; break; }
			if (id == -1) return {false, matrix()};
			if (id != i) {
				swaps.push_back({id, i});
				for (int j = 0; j < n; j++) std::swap(a[i][j], a[id][j]);
			}
			a[i][i] =  T(1) / a[i][i];
			for (int j = 0; j < n; j++) if (j != i) a[i][j] *= a[i][i];
			for (int j = 0; j < n; j++) if (j != i) {
				for (int k = 0; k < n; k++) if (k != i) a[j][k] -= a[j][i] * a[i][k];
				a[j][i] *= -a[i][i];
			}
		}
		for (int i = swaps.size(); i--; ) 
			for (int j = 0; j < n; j++) std::swap(a[j][swaps[i].first], a[j][swaps[i].second]);
		return {true, a};
	}
	T det() const {
		assert(n == m);
		matrix a = *this;
		T res = 1;
		for (int i = 0; i < n; i++) {
			int id = -1;
			for (int j = i; j < n; j++) if (a[j][i] != T(0)) { id = j; break; }
			if (id == -1) return 0;
			if (id != i) {
				res = -res;
				for (int j = i; j < n; j++) std::swap(a[id][j], a[i][j]);
			}
			res *= a[i][i];
			
			T t = T(1) / a[i][i];
			for (int j = i; j < n; j++) a[i][j] *= t;
			for (int j = i + 1; j < n; j++) {
				auto s = a[j][i];
				for (int k = i; k < n; k++) a[j][k] -= a[i][k] * s;
			}
		}
		return res;
	}
};
