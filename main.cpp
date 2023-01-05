#include "bits/stdc++.h"
using namespace std;

#include "debug.h"

const int inf = 1e9;

std::mt19937 rng(1);

double randDouble(double lo = 0.0, double hi = 1.0) {
    std::uniform_real_distribution<> dis(lo, hi);
    return dis(rng);
}


struct Perferences{
	int n; // number of pairs of men and women

	vector<vector<int>> man, woman, mn, wmn;
	// man[i][j] -> ith man to jth woman preference score
	// high score means high perference, -1 means this marriage can't possible..
	// mn[i] -> Perferences list of ith man

	Perferences(int _n) : n(_n), man(_n, vector<int>(_n)), woman(_n, vector<int>(_n)), mn(_n), wmn(_n) {}

	void random_gen_SMP(){
		// Random generator for Stable Marriage Problem
		for(auto &v : man){
			iota(v.begin(), v.end(), 0);
			shuffle(v.begin(), v.end(), rng);
		}

		for(auto &v : woman){
			iota(v.begin(), v.end(), 0);
			shuffle(v.begin(), v.end(), rng);
		}

		gen_mn_wmn();
	}

	void random_gen_SMPT(){
		// Random generator for Stable Marriage Problem with ties

		for(auto &v : man){
			for(auto &x : v) x = rng()%n;

			auto z = v;
			sort(z.begin(), z.end());
			z.erase(unique(z.begin(), z.end()), z.end());

			for(auto &x : v) x = lower_bound(z.begin(), z.end(), x) - z.begin();
		}

		for(auto &v : woman){
			for(auto &x : v) x = rng()%n;

			auto z = v;
			sort(z.begin(), z.end());
			z.erase(unique(z.begin(), z.end()), z.end());

			for(auto &x : v) x = lower_bound(z.begin(), z.end(), x) - z.begin();
		}

		gen_mn_wmn();
	}

	void make_incomplete(double p){
		for(auto &v : man){
			for(auto &x : v) if(randDouble() < p) x = -1;
		}

		for(auto &v : woman){
			for(auto &x : v) if(randDouble() < p) x = -1;
		}
	}


	void random_gen_SMPTI(double p = 0.2){ // low p means low incomplete pairs
		// Random generator for Stable Marriage Problem with ties with Incomplete List
		// woman[i][j] = -1 means ith woman can't prefer jth man

		random_gen_SMPT();
		make_incomplete(p);
		gen_mn_wmn();
	}

	void random_gen_SMI(double p = 0.2){
		// Random generator for Stable Marriage Problem with Incomplete List
		random_gen_SMP();
		make_incomplete(p);
		gen_mn_wmn();
	}

	void gen_mn_wmn(){
		for(int i = 0; i < n; ++i){
			auto &v = mn[i], &u = man[i];
			for(int j = 0; j < n; ++j) if(u[j] != -1){
				v.emplace_back(j);
			}
			sort(v.begin(), v.end(), [&](int a, int b){
				return u[a] > u[b];
			});
		}

		for(int i = 0; i < n; ++i){
			auto &v = wmn[i], &u = woman[i];
			for(int j = 0; j < n; ++j) if(u[j] != -1){
				v.emplace_back(j);
			}
			sort(v.begin(), v.end(), [&u](int a, int b){
				return u[a] > u[b];
			});
		}
	}

};

struct Grader{

	int n;
	const Perferences &p;
	// Matching m;
	const vector<vector<int>> &mp, &wp;

	// Grader(const Perferences &_p, Matching &_m) : n(_p.n), p(_p), m(_m), mp(p.man), wp(p.woman) {}
	Grader(const Perferences &_p) : n(_p.n), p(_p), mp(p.man), wp(p.woman) {}

	double egalitarian_score(vector<int> &match){
		int sc = 0;

		for(int m = 0, w; m < n; ++m){
			w = match[m];
			sc += mp[m][w] - wp[w][m];
		}

		return ((double)abs(sc))/n;
	}

	int stable_pairs(vector<int> &match){
		int sc = 0;
		vector<bool> g(n, 1);

		// dbg(mp, wp)

		for(int m1 = 0; m1 < n; ++m1){
			for(int m2 = m1 + 1; m2 < n; ++m2){
				int w1 = match[m1], w2 = match[m2];

				if(mp[m1][w1] < mp[m1][w2] && wp[w2][m2] < wp[w2][m1]){
					// m1, w2 will marry
					// dbg(m1, m2, w1, w2)
					g[m1] = g[m2] = 0;
					break;
				}
			}
		}

		return accumulate(g.begin(), g.end(), (int)0);
	}
};


struct Matching{

	int n;
	const Perferences &p;
	vector<int> v, vw;
	const vector<vector<int>> &mp, &wp;

	Grader G;

	Matching(const Perferences &_p) : n(_p.n), p(_p), v(n, -1), vw(n, -1), mp(_p.man), wp(_p.woman), G(p) {}

	double GradeMe(){
		double egalitarian_score = G.egalitarian_score(v);
		int stable_pairs = G.stable_pairs(v);

		dbg(egalitarian_score, stable_pairs);
		return - egalitarian_score + stable_pairs;
	}

};

struct GSA : Matching {
	const vector<vector<int>> &mn, &wmn;

	GSA(const Perferences &p) : Matching(p), mn(p.mn), wmn(p.wmn) {
		vector<int> cm(n); // current mapping of man in terms of index

		for(int _m = 0, m = 0, w; m < n; ++_m, m = _m){

			while(m != -1 && cm[m] < mn[m].size()){
				w = mn[m][cm[m]++];

				if(vw[w] == -1 || wp[w][vw[w]] < wp[w][m]){
					v[m] = w;
					swap(m, vw[w]);
				}
			}
		}
		GradeMe();
	}
};

struct HILL : Matching {
	double mxscore;
	
	HILL(const Perferences &p) : Matching(p), mxscore(0) {
		// anneling
		for(int i = 0; i < 1000; ++i){
			auto cv = random_stable();
			double sc = FitnessScore(cv);
			if(sc > mxscore){
				mxscore = sc;
				swap(v, cv);
			}
		}

		GradeMe();
	}

	double FitnessScore(vector<int> &p){
		return (double)G.stable_pairs(p) - G.egalitarian_score(p);
	}


	// finding random stable matching
	// avg(O(n^2))
	vector<int> random_stable(){

		vector<int> v(n), w;
		iota(v.begin(), v.end(), 0);
		w = v;
		shuffle(v.begin(), v.end(), rng);
		shuffle(w.begin(), w.end(), rng);

		while(w.size()){
			for(int i = 0; i < w.size(); ++i){
				int m1 = w[i], w1 = v[m1];
				bool ok = 0;
				for(int j = i + 1; j < w.size(); ++j){
					int m2 = w[j], w2 = v[m2];
					
					if(mp[m1][w1] < mp[m1][w2] && wp[w2][m2] < wp[w2][m1]){
						ok = 1;
						swap(v[m1], v[m2]);
						break;
					}
				}
				if(!ok){
					swap(w.back(), w[i]);
					w.pop_back(); --i;
				}
			}
		}

		// vector<vector<int>> u(n);
		// for(int m1 = 0; m1 < n; ++m1){
		// 	for(int m2 = m1 + 1; m2 < n; ++m2){
		// 		int w1 = v[m1], w2 = v[m2];

		// 		if(mp[m1][w1] < mp[m1][w2] && wp[w2][m2] < wp[w2][m1]){
		// 			// m1, w2 will marry
		// 			u[m1].emplace_back(m2);
		// 			u[m2].emplace_back(m1);
		// 		}
		// 	}
		// }
		// dbg(u)
		// vector<int> z;
		// for(int i = 0; i < n; ++i) if(u[i].size()) z.emplace_back(i);

		// int rep = n*n;


		// for(int _ = 0; _ < rep; ++_){
		// 	dbg(z, v)
		// 	int _q = rng()%z.size();
		// 	int m1 = z[_q];

		// 	auto &x = u[m1];
		// 	dbg(m1, x)
		// 	if(x.size() == 0){
		// 		z.erase(z.begin() + _q); continue;
		// 	}

		// 	int q = rng()%x.size();
		// 	int m2 = x[q], w1 = v[m1], w2 = v[m2];
		// 	dbg(q, m2, w1, w2)

		// 	if(mp[m1][w1] > mp[m1][w2] || wp[w2][m2] > wp[w2][m1]){
		// 		x.erase(x.begin() + q);
		// 		dbg(x)
		// 		continue;
		// 	}

		// 	x.resize(q);

		// 	auto &y = u[m2];
		// 	while(y.back() != m1) y.pop_back();
		// 	y.pop_back();

		// 	swap(v[m1], v[m2]);
		// }

		return v;
	}
};

struct Genetic : Matching {

	int Popfactor, PopulationSize;
	// Population Size = Popfactor * n;


	using Gen = vector<vector<int>>; // Generation
	Gen P; // Population
	const int MR1; // Mutation Rate 1
	double mxscore;


	Genetic(const Perferences &p, int Popfactor = 20) : Matching(p), Popfactor(Popfactor), PopulationSize(Popfactor*n), P(PopulationSize, vector<int>(n)), MR1((n + 19)/20) {
		mxscore = 0;

		int rep = 100;
		set<vector<int>> Z;
		for(int _ = 0; _ < rep; ++_){
			P = Gen(PopulationSize, vector<int>(n));
			Init();

			for(auto &x : Z) P.emplace_back(x);

			for(int _z = 0; _z < 30; ++_z){
				Gen Childs = CrossOver(P);
				for(auto &x : Childs){
					if(rng()&1) Mutation1(x);
					else{
						for(int i = 0; i < n; ++i) if(!Mutation2(x)) break;
					}
					P.emplace_back(x);
				}

				double F = 0;
				vector<double> cs;

				for(auto &x : P){
					auto z = FitnessScore(x);
					F += z;
					if(z > mxscore){
						mxscore = z;
						v = x;
					}
					cs.emplace_back(z*PopulationSize);
				}


				Gen nP; // new Population
				for(int i = 0; i < P.size(); ++i){
					if(randDouble() < cs[i]/F){
						nP.emplace_back(P[i]);		
					}
				}

				swap(P, nP);
				shuffle(P.begin(), P.end(), rng);

			}

			vector<pair<double, int>> vp;
			for(int i = 0; i < P.size(); ++i){
				vp.emplace_back(FitnessScore(P[i]), i);
			}

			sort(vp.rbegin(), vp.rend());
			for(int i = 0, zz = min(PopulationSize/10, (int)vp.size()); i < zz; ++i){
				Z.emplace(P[vp[i].second]);
			}

			GradeMe();
		}

	}

	void Init(){
		for(auto &x : P){
			iota(x.begin(), x.end(), 0);
			shuffle(x.begin(), x.end(), rng);
		}
	}

	Gen CrossOver(const Gen &P){

		int N = P.size()/2*2;
		Gen K(N, vector<int>(n, -1)); // Kids

		for(int i = 0; i < N; i += 2){
			auto &p1 = P[i], &p2 = P[i + 1];
			vector<int> p1i(n);
			for(int i = 0; i < n; ++i){
				p1i[p1[i]] = i; 
			}
			auto &k1 = K[i], &k2 = K[i + 1];

			int ptr = rng()%n;

			while(k1[ptr] == -1){
				k1[ptr] = p1[ptr];
				k2[ptr] = p2[ptr];
				ptr = p1i[p2[ptr]];
			}
			// dbg(ptr, k1, k2, p1, p2);

			for(int j = 0; j < n; ++j) if(k1[j] == -1){
				k1[j] = p2[j];
				k2[j] = p1[j];
			}
			// dbg(ptr, k1, k2, p1, p2);
		}


		return K;
	}

	void Mutation1(vector<int> &child){
		for(int _ = 0; _ < MR1; ++_){
			int a = rng()%n;
			int b = rng()%n;
			swap(child[a], child[b]);
		}
	}

	bool Mutation2(vector<int> &child){
		// find a random unstable pair and swap them
		vector<int> mans(n);
		iota(mans.begin(), mans.end(), 0);
		shuffle(mans.begin(), mans.end(), rng);

		for (std::vector<int>::iterator i = mans.begin(); i != mans.end(); ++i){
			int m1 = *i, w1 = child[m1];
			for (std::vector<int>::iterator j = i + 1; j != mans.end(); ++j){
				int m2 = *j, w2 = child[m2];
				if(mp[m1][w1] < mp[m1][w2] && wp[w2][m2] < wp[w2][m1]){
					swap(child[m1], child[m2]);
					return 1;
				}
			}
		}

		return 0;
	}

	double FitnessScore(vector<int> &p){
		return (double)G.stable_pairs(p) - G.egalitarian_score(p);
	}

};



signed main(int argc, char const *argv[]){
	if(argc > 1){
		rng.seed(stoi(argv[1]));
		srand(stoi(argv[1]));
	}else{
		srand(1);
	}

	Perferences P(50);
	P.random_gen_SMP();

	GSA G(P);
	Genetic G1(P);
	HILL G2(P);
}